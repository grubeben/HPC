/* (C) Jesper Larsson Traff, October 2022 */
/* Alltoall algorithms for fully connected networks */
/* Example code for HPC 2022, see script Section 7.2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <mpi.h>

// #include "tuw_alltoall.h"

/* Algorithm selection done at compile time */
#ifndef ALLTOALL
#define ALLTOALL ReduceBcast
// #define ALLTOALL MPI_Allreduce
#endif

#define ALLTOALLTAG 7777

// Benchmarking parameters
#define WARMUP 8
#define REPEAT 41     // odd number to easier compute the median
#define MEDIAN_IDX 20 // index of 21st element in runtime array
#define Z_VALUE 2.57  // corresponds to 99% confidence interval (1.64...90%, 1.96...95%)
#define MICRO 1000000.0

#define TUW_TYPE MPI_FLOAT
typedef float tuwtype_t;

#define FLIP(_fac, _mul1, _mul2) (_fac = ((_fac == _mul1) ? (_mul2) : (_mul1)))

int ReduceBcast(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm)
{
    // int rank, size;

    // is it sendcount or recvcount ?
    MPI_Reduce(sendbuf, recvbuf, sendcount, sendtype, MPI_MAX, 0, comm);

    // Barrier necessary ? ----> No,
    // all collective communications are non-blocking
    //MPI_Barrier(comm);
    //MPI_Barrier(comm);

    // warum macht er (char*) bei sendbuf ?
    // sendcount als reference ?
    MPI_Bcast(recvbuf, sendcount, sendtype, 0, comm);

    return MPI_SUCCESS;
}

void statistics(double *runtime, int c, int size, char x)
{

    int i, j;
    double tuwavg, tuwmed, tuwmin, tmp;

    // sorting of the runtimes
    tuwavg = 0;
    for (i = 0; i < REPEAT; ++i)
    {
        tuwavg += runtime[i];
        for (j = i + 1; j < REPEAT; ++j)
        {

            if (runtime[i] > runtime[j])
            {

                tmp = runtime[i];
                runtime[i] = runtime[j];
                runtime[j] = tmp;
            }
        }
    }

    tuwavg /= REPEAT;
    tuwmin = runtime[0];
    tuwmed = runtime[MEDIAN_IDX];

    double stddev = 0;
    for (i = 0; i < REPEAT; ++i)
    {
        stddev += (runtime[i] - tuwavg) * (runtime[i] - tuwavg);
    }
    stddev = sqrt(stddev / REPEAT);

    // x(+/-)t*(s/√n)

    double lower, upper;
    upper = tuwavg + Z_VALUE * (stddev / sqrt(REPEAT));
    lower = tuwavg - Z_VALUE * (stddev / sqrt(REPEAT));

    printf("%c %8d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
           x, c, tuwavg * MICRO, stddev * MICRO, tuwmed * MICRO, tuwmin * MICRO, lower * MICRO, upper * MICRO);
}

int main(int argc, char *argv[])
{
    int rank, size;
    int count, c;

    tuwtype_t *sendbuf, *recvbuf, *testbuf;

    int i;
    int r, t;

    double start, stop;
    double runtime[REPEAT], runtimeMPI[REPEAT];

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of elements
    count = 1000000;

    // allocate memory
    sendbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf != NULL);
    recvbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf != NULL);
    testbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf != NULL);

    // initialize data (fill arrays)
    for (i = 0; i < count; i++)
        sendbuf[i] = (tuwtype_t)i;
    for (i = 0; i < count; i++)
        recvbuf[i] = (tuwtype_t)0;
    for (i = 0; i < count; i++)
        testbuf[i] = (tuwtype_t)0;

    // Manual "correctness test"
    // if (rank == 0)
    // {
    //     printf("\nReduceBcast:\n");
    //     for (int i = 0; i < count; ++i)
    //     {
    //         printf("%.2f\n", recvbuf[i]);
    //     }
    //     printf("\nMPI_Allreduce:\n");
    //     for (int i = 0; i < count; ++i)
    //     {
    //         printf("%.2f\n", testbuf[i]);
    //     }
    // }

    // "correctness test": compare against result from library function
    ReduceBcast(sendbuf, count, TUW_TYPE, recvbuf, count, TUW_TYPE, MPI_COMM_WORLD);
    MPI_Allreduce(sendbuf, testbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    for (i = 0; i < count; i++)
        assert(recvbuf[i] == testbuf[i]);

    //   if (rank==0) {
    //     fprintf(stderr,"Alltoall, size=%d\\\\\n",size);
    //     fprintf(stderr,"count & m (count) & m (Bytes) & avg & min \\\\\n");
    //   }
    if (rank==0)
    printf("%s %6s %8s %8s %8s %8s %8s %8s\n", "ver", "N", "avg", "std", "med", "min", "lower", "upper");

    int f = 10;
    // c*=f,FLIP(f,2,5)
    // --> do the c*=f and then execute FLIP(f,2,5),
    // which writes into the f the result FLIP(...)
    // FLIP(_fac, a, b) (_fac = (_fac==a) ? b : a )
    // FLIP(f,2,5) increases c by powers of 2 and powers of 5 alternatingly (c*=2, c*=5, c*=2, c*=5, ...)
    for (int c = 1; c <= count; c *= f, FLIP(f, 2, 10))
    {
        // if (c*size>size*count) break;

        for (r = 0, t = 0; r < WARMUP + REPEAT; r++)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD); // good practice -> barrier does not work super reliably

            start = MPI_Wtime();
            ReduceBcast(sendbuf, c, TUW_TYPE, recvbuf, c, TUW_TYPE, MPI_COMM_WORLD);
            stop = MPI_Wtime();

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            if (r >= WARMUP)
            {
                runtime[t] = stop - start;
            }

            start = MPI_Wtime();
            MPI_Allreduce(sendbuf, testbuf, c, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
            stop = MPI_Wtime();

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            if (r < WARMUP)
                continue;
            runtimeMPI[t++] = stop - start;
        }
        // reduce the timings of each processor using always the max value
        MPI_Allreduce(MPI_IN_PLACE, runtime, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);    //
        MPI_Allreduce(MPI_IN_PLACE, runtimeMPI, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); //
                                                                                         // REICHT HIER EIN NORMALES REDUCE ?

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // proc 0 computes statistics of runtimes
        if (rank == 0)
        {
            statistics(runtime, c, size, 'c');    // 'c' for custom
            statistics(runtimeMPI, c, size, 'r'); // 'r' for reference
        }
    }

    /* -------- TO DO --------
    - run on hydra (with all process configurations (nodes, count) given on sheet)
    - explain benchmark setup
    - Plots überdenken

      ---------- Q&A ---------
    - passt die Sache mit dem Confidence interval ?
    - Hydra erlaubt nicht mehr als ?? processes (error bei 16 und 32)
    - Adrian's code um Faktor 2-4 schneller
    - mit FLOATs langsamer als mit DOUBLEs !?
    - timings werden länger mit steigender prozessoranzahl wtf!? DONE
            (--> sollte ja eh so sein)
    - first run of MPI_Allreduce produces ~ 10x larger results DONE 
            (--> wir haben 'continue' falsch verstanden)
    */

    MPI_Finalize();

    free(sendbuf);
    free(recvbuf);
    free(testbuf);

    return 0;
}
