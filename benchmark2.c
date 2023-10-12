#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#include "allreduce.h"

// Benchmarking parameters
#define WARMUP 2
#define REPEAT 21    // odd number to easier compute the median
#define Z_VALUE 1.96  // corresponds to 99% confidence interval (1.64...90%, 1.96...95%, 2.57...99%)
#define MICRO 1000000.0

#define FLIP(_fac, _mul1, _mul2) (_fac = ((_fac == _mul1) ? (_mul2) : (_mul1)))
#define MIN(X, Y) (((X) < (Y)) ? (x) : (Y))

void statistics(double *runtime, int c, int blocksize, int size, int repeat, int median_idx, char x)
{
    int i, j;
    double tuwavg, tuwmed, tuwmin, tmp;

    // sorting of the runtimes
    tuwavg = 0;
    for (i = 0; i < repeat; ++i)
    {
        tuwavg += runtime[i];
        for (j = i + 1; j < repeat; ++j)
        {
            if (runtime[i] > runtime[j])
            {
                tmp = runtime[i];
                runtime[i] = runtime[j];
                runtime[j] = tmp;
            }
        }
    }

    tuwavg /= repeat;
    tuwmin = runtime[0];
    tuwmed = runtime[median_idx];

    double stddev = 0;
    for (i = 0; i < repeat; ++i)
    {
        stddev += (runtime[i] - tuwavg) * (runtime[i] - tuwavg);
    }
    stddev = sqrt(stddev / repeat);

    double lower, upper;
    upper = tuwavg + Z_VALUE * (stddev / sqrt(repeat));
    lower = tuwavg - Z_VALUE * (stddev / sqrt(repeat));

    printf("%c %d %d %d %.4f %.4f %.4f %.4f %.4f %.4f\n",
           x, size, c, blocksize, tuwavg * MICRO, stddev * MICRO, tuwmed * MICRO, tuwmin * MICRO, lower * MICRO, upper * MICRO);
}

int counts[10] = {1, 2, 100, pow(2, 9), pow(2, 12), 10000, pow(2, 16), pow(2, 19), 1000000, 10000000};

int main(int argc, char *argv[])
{
    int repeat = REPEAT;
    int median_idx;
    int count = 10000000; // number of elements
    int rank, size;
    int c; // running number of elements
    int i; // loop index
    int r, t;
    int testBlockSize = 3 * count / 4; // blocksize for assert

    tuwtype_t *sendbuf, *sendbuf1, *testbuf, *recvbuf1, *recvbuf2, *recvbuf3, *recvbuf4, *recvbuf5;

    double start, start1, start2, start3, start4, start5, stop, stop1, stop2, stop3, stop4, stop5;
    double runtime1[REPEAT], runtime2[REPEAT], runtime3[REPEAT], runtime4[REPEAT], runtime5[REPEAT], runtimeMPI[REPEAT];

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // allocate memory
    sendbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf != NULL);
    sendbuf1 = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf1 != NULL);
    testbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf != NULL);
    recvbuf1 = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf1 != NULL);
    recvbuf2 = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf2 != NULL);
    recvbuf3 = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf3 != NULL);
    recvbuf4 = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf4 != NULL);
    recvbuf5 = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf5 != NULL);

    // initialize data (fill arrays)
    for (i = 0; i < count; i++)
    {
        sendbuf[i] = (tuwtype_t)rank;
        sendbuf1[i] = (tuwtype_t)rank;
        testbuf[i] = (tuwtype_t)0;
        recvbuf1[i] = (tuwtype_t)0;
        recvbuf2[i] = (tuwtype_t)rank;
        recvbuf3[i] = (tuwtype_t)rank;
        recvbuf4[i] = (tuwtype_t)rank;
        recvbuf5[i] = (tuwtype_t)rank;
    }

    // COORECTNESS TEST EXERCISE 1-5
    allReduce1(sendbuf1, count, recvbuf1, count, MPI_COMM_WORLD);
    allReduce2(recvbuf2, count, testBlockSize, size, rank);
    allReduce3(recvbuf3, count, testBlockSize, size, rank);
    allReduce4(recvbuf4, count, testBlockSize, size, rank);
    allReduce5(recvbuf5, count, testBlockSize, size, rank);

    // MPI library function reference
    MPI_Allreduce(sendbuf, testbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    for (i = 0; i < count; i++)
    {
        assert(recvbuf1[i] == testbuf[i]);
        assert(recvbuf2[i] == testbuf[i]);
        assert(recvbuf3[i] == testbuf[i]);
        assert(recvbuf4[i] == testbuf[i]);
        assert(recvbuf5[i] == testbuf[i]);
    }

    // Initiate header for output
    if (rank == 0)
        printf("%s %s %s %s %s %s %s %s %s %s\n", "ver", "size", "count", "b", "avg", "std", "med", "min", "lower", "upper");

    // Control loop settings
    //int fc = 10;
    int fb = 10;

    //for (int c = 1; c < count; c *= fc, FLIP(fc, 2, 10))
    for (int i = 0; i < 10; ++i)
    {
        c = counts[i];
        repeat -= 2;
        median_idx = ((int)(repeat/2));

        /* blockSize = 10 100 1000 10000 100000 1000000 10000000*/
        for (int blockSize = 10; blockSize <= c; blockSize *= 10)
        {
            if (blockSize > c) 
                blockSize = c;

            /*
            if (blockSize + c/4+1 > c)
                blockSize = c;
            */

            for (r = 0, t = 0; r < WARMUP + repeat; r++)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD); // good practice -> barrier does not work super reliably

                ///////////////////////////////////
                // Insert desired function below //
                ///////////////////////////////////
                start = MPI_Wtime();
                MPI_Allreduce(sendbuf, testbuf, c, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
                stop = MPI_Wtime();
                

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);

                start1 = MPI_Wtime();
                //allReduce1(sendbuf1, c, recvbuf1, c, MPI_COMM_WORLD);
                stop1 = MPI_Wtime();
                

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);

                start2 = MPI_Wtime();
                allReduce2(recvbuf2, c, blockSize, size, rank);
                stop2 = MPI_Wtime();

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);

                start3 = MPI_Wtime();
                //allReduce3(recvbuf3, c, blockSize, size, rank);
                stop3 = MPI_Wtime();

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);

                start4 = MPI_Wtime();
                //allReduce4(recvbuf4, c, blockSize, size, rank);
                stop4 = MPI_Wtime();

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);

                start5 = MPI_Wtime();
                //allReduce5(recvbuf5, c, blockSize, size, rank);
                stop5 = MPI_Wtime();

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);

                ////////////////////////////

                if (r >= WARMUP)
                {
                    runtime1[t] = stop1 - start1;
                    runtime2[t] = stop2 - start2;
                    runtime3[t] = stop3 - start3;
                    runtime4[t] = stop4 - start4;
                    runtime5[t] = stop5 - start5;
                }

                if (r < WARMUP)
                    continue;  
                runtimeMPI[t++] = stop - start;              
            }
            // reduce the timings of each processor using always the max value
            MPI_Allreduce(MPI_IN_PLACE, runtime1, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, runtime2, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, runtime3, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, runtime4, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, runtime5, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, runtimeMPI, t, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            // node 0 computes statistics of runtimes
            if (rank == 0)
            {
                //statistics(runtime1, c, blockSize, size, repeat, median_idx, '1'); // 'c' for custom
                statistics(runtime2, c, blockSize, size, repeat, median_idx, '2'); // 'c' for custom
                //statistics(runtime3, c, blockSize, size, repeat, median_idx, '3'); // 'c' for custom
                //statistics(runtime4, c, blockSize, size, repeat, median_idx, '4'); // 'c' for custom
                //statistics(runtime5, c, blockSize, size, repeat, median_idx, '5'); // 'c' for custom
                statistics(runtimeMPI, c, blockSize, size, repeat, median_idx, 'r'); // 'r' for reference
            }
        }
    }

    MPI_Finalize();

    free(sendbuf);
    free(sendbuf1);
    free(testbuf);
    free(recvbuf1);
    free(recvbuf2);
    free(recvbuf3);
    free(recvbuf4);
    free(recvbuf5);

    return 0;
}