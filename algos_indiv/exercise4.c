#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include <assert.h>

#include <mpi.h>

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

int treeduce(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int childL = 2 * rank + 1;
    int childR = 2 * rank + 2;
    int childLexists = (childL < size);
    int childRexists = (childR < size);

    double d = ceil(log2(size + 1) - 1); // binary tree depth
    int perfect_size = 0;
    for (int i = 0; i <= (int)d; ++i)
    {
        perfect_size += pow(2, i);
    }

    int parent;
    // determine parent based on even/oddness of rank
    if (rank % 2 == 0)
    {
        parent = (rank - 2) / 2;
    }
    else
    {
        parent = (rank - 1) / 2;
    }

    // printf("PROCESS %d with parent: %d and depth %f | childL: %d | childR: %d\n", rank, parent, d, childL, childR);

    int trunc_chunk_idx = floor(count / blockSize) * blockSize; // index of first truncated chunk
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbuf_trunc = (tuwtype_t *)malloc(count % blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;
    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            // printf("ROOTNODE %d in line 104\n", rank);
            if (childLexists)
            {
                MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, MPI_DOUBLE, MPI_MAX);
            }
            // MPI_Barrier(MPI_COMM_WORLD);
            if (childRexists)
            {
                MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, MPI_DOUBLE, MPI_MAX);
            }
            // printf("ROOTNODE %d in line 110 \n", rank);
        }
        if (count % blockSize != 0)
        {
            // printf("ROOTNODE %d in line 114 leftover block\n", rank);
            if (childLexists)
            {
                MPI_Recv(tmpbuf_trunc, count % blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, MPI_MAX);
            }
            // MPI_Barrier(MPI_COMM_WORLD);
            if (childRexists)
            {
                MPI_Recv(tmpbuf_trunc, count % blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, MPI_MAX);
            }
        }
    }
    else if ((!childRexists) && (!childLexists)) // check if this is a leaf-node
    {
        // printf("send only node start\n");
        // printf("LEAF %d in line 125\n", rank);
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(sendbuf + b, blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD);
        }
        if (count % blockSize != 0)
        {
            MPI_Send(sendbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD);
        }
        // printf("LEAF %d in line 134\n", rank);
    }
    else
    {
        if (count >= blockSize)
        {
            // printf("MIDDELNODE %d in line 153\n", rank);
            for (int b = 0; b < trunc_chunk_idx; b += blockSize)
            {
                // printf("MIDDELNODE %d in line 156\n", rank);
                if (childLexists)
                {
                    // printf("MIDDELNODE %d in line 159\n", rank);
                    MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, MPI_DOUBLE, MPI_MAX);
                }
                if (childRexists)
                {
                    // printf("MIDDELNODE %d in line 165 \n", rank);
                    MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, MPI_DOUBLE, MPI_MAX);
                }
                // printf("MIDDELNODE %d in line 169 \n", rank);
                MPI_Send(recvbuf + b, blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD);
                // printf("MIDDELNODE %d in line 171 \n", rank);
            }
            // printf("MIDDELNODE %d in line 176 \n", rank);
        }
        if (count % blockSize != 0)
        {
            // printf("MIDDELNODE %d in line 180 leftover block \n", rank);
            if (childLexists)
            {
                MPI_Recv(tmpbuf_trunc, count % blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, MPI_MAX);
            }
            // MPI_Barrier(MPI_COMM_WORLD);
            if (childRexists)
            {
                MPI_Recv(tmpbuf_trunc, count % blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, MPI_MAX);
            }
            // MPI_Barrier(MPI_COMM_WORLD);
            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD);
            // printf("MIDDELNODE %d in line 188 leftover block done \n", rank);
        }
    }
    return MPI_SUCCESS;
}

// binary tree broadcast
int treecast(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int childL = 2 * rank + 1;
    int childR = 2 * rank + 2;
    int childLexists = (childL < size);
    int childRexists = (childR < size);

    double d = ceil(log2(size + 1) - 1); // binary tree depth
    int perfect_size = 0;
    for (int i = 0; i <= (int)d; ++i)
    {
        perfect_size += pow(2, i);
    }

    int parent;
    // determine parent based on even/oddness of rank
    if (rank % 2 == 0)
    {
        parent = (rank - 2) / 2;
    }
    else
    {
        parent = (rank - 1) / 2;
    }

    // printf("PROCESS %d with parent: %d and depth %f | childL: %d | childR: %d\n", rank, parent, d, childL, childR);

    int trunc_chunk_idx = floor(count / blockSize) * blockSize;

    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            if (childLexists)
                MPI_Send(recvbuf + b, blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD);
            if (childRexists)
                MPI_Send(recvbuf + b, blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD);
        }

        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            if (childLexists)
                MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD);
            if (childRexists)
                MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD);
        }
    }

    else if ((!childRexists) && (!childLexists)) // check if this is a leaf-node
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(recvbuf + b, blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Recv(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        if (count >= blockSize)
        {
            // printf("%d\n",rank);
            for (int b = 0; b < trunc_chunk_idx; b += blockSize)
            {
                // printf("procc %d loop %d\n", rank, b);
                //  MPI_Sendrecv(sendbuf + b - blockSize, blockSize, MPI_DOUBLE, rank + 1, 0, recvbuf + b, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(recvbuf + b, blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (childLexists)
                    MPI_Send(recvbuf + b, blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD);
                if (childRexists)
                {
                    MPI_Send(recvbuf + b, blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD);
                }
            }
        }
        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Recv(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (childLexists)
                MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD);
            if (childRexists)
                MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD);
        }
    }

    return MPI_SUCCESS;
}

int PipelinedTreeduceTreecast(void *sendbuf, void *recvbuf, int count, int blockSize, int size, int rank)
{
    if (size == 1)
    {
        memcpy(recvbuf, sendbuf, count * sizeof(double));
        // could be done more elegantly with pointers.
        // Also the datatype should not be hardcoded
        return MPI_SUCCESS;
    }

    treeduce(sendbuf, recvbuf, count, blockSize, size, rank);
    treecast(sendbuf, recvbuf, count, blockSize, size, rank);

    return MPI_SUCCESS;
}

int main(int argc, char *argv[])
{
    int rank, size, count, blockSize = 3, i;
    tuwtype_t *sendbuf, *recvbuf, *testbuf;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    count = 10;

    sendbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf != NULL);
    recvbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf != NULL);
    testbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf != NULL);

    srand(time(NULL));

    for (i = 0; i < count; i++)
        sendbuf[i] = (tuwtype_t)rand();
    for (i = 0; i < count; i++)
        recvbuf[i] = (tuwtype_t)-1;
    for (i = 0; i < count; i++)
        testbuf[i] = (tuwtype_t)-1;

    // "correctness test": compare against result from library function

    /* do we have to pass size and rank, no smooth way to identify within method? */
    // reduce(sendbuf, testbuf, count, size, rank);
    // bcast(sendbuf, testbuf, count, size, rank);
    PipelinedTreeduceTreecast(sendbuf, testbuf, count, blockSize, size, rank);
    // treecast(sendbuf, testbuf, count, blockSize, size, rank);
    // treeduce(sendbuf, testbuf, count, blockSize, size, rank);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("reference  custom\n");
        for (i = 0; i < count; i++)
        {
            // assert(recvbuf[i] == testbuf[i]);
            printf("%f, %f\n", recvbuf[i], testbuf[i]);
        }
    }

    MPI_Finalize();

    free(sendbuf);
    free(recvbuf);
    free(testbuf);

    return 0;
}
