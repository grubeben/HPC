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
