#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include <assert.h>

#include <mpi.h>

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int doublyPipelinedAlltreeduce(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int childL = 2 * rank + 1;
    int childR = 2 * rank + 2;
    int childLexists = (childL < size);
    int childRexists = (childR < size);

    int d_fac = ceil(log2(size + 1) - 1); // binary tree depth
    int d = d_fac * blockSize;
    int di; // node depth *BlockSize

    // find node depth
    double size_run = 0;
    for (int d_run = 0; d_run < d_fac + 1; d_run++)
    {
        size_run += pow(2, d_run);
        if (rank < size_run)
        {
            di = d_run * blockSize;
            break;
        }
    }
    // printf("Node %d is on depth %f \n",rank,di );

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

    int last_block_idx = floor(count / blockSize) * blockSize; // index of first truncated chunk
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbuf_trunc = (tuwtype_t *)malloc(count % blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;

    if (rank == 0) // root has no phase1
    {
        // copy sendbuf into recvbuf
        MPI_Reduce_local(sendbuf, recvbuf, count, MPI_DOUBLE, MPI_MAX);
        for (int block_idx = 0; block_idx < last_block_idx + di + d + blockSize; block_idx += blockSize)
        {
            // which phase is root in?
            if (block_idx < d && block_idx < last_block_idx) // phase 1 - until first non-zero block arrives at root
            {
                if (childLexists)
                    MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, childL, 0, recvbuf, 0, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (childRexists)
                    MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, childR, 0, recvbuf, 0, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else if (block_idx = d) // induce offset
            {
                if (childLexists)
                {
                    MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, childL, 0, tmpbuf, blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + (block_idx - d), blockSize, MPI_DOUBLE, MPI_MAX);
                }
                if (childRexists)
                {
                    MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, childR, 0, tmpbuf, blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + (block_idx - d), blockSize, MPI_DOUBLE, MPI_MAX);
                }
            }
            else if (block_idx > d && block_idx <= last_block_idx + d) // phase 2 - from receiving second block to receiving last non-zero blocks;
            {
                if (childLexists)
                {
                    MPI_Sendrecv(recvbuf + (block_idx - d - blockSize), blockSize, MPI_DOUBLE, childL, 0, tmpbuf, blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + (block_idx - d), blockSize, MPI_DOUBLE, MPI_MAX);
                }
                if (childRexists)
                {
                    MPI_Sendrecv(recvbuf + (block_idx - d - blockSize), blockSize, MPI_DOUBLE, childR, 0, tmpbuf, blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + (block_idx - d), blockSize, MPI_DOUBLE, MPI_MAX);
                }
            }
            else if (block_idx = last_block_idx + d + blockSize) // redeem offset
            {
                if (childLexists)
                {
                    MPI_Sendrecv(recvbuf + (block_idx - d - blockSize), blockSize, MPI_DOUBLE, childL, 0, recvbuf, 0, MPI_DOUBLE, childL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                if (childRexists)
                {
                    MPI_Sendrecv(recvbuf + (block_idx - d - blockSize), blockSize, MPI_DOUBLE, childR, 0, recvbuf, 0, MPI_DOUBLE, childR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    else if ((!childRexists) && (!childLexists)) // check if this is a leaf-node
    {
        // printf("send only node start\n");
        // printf("LEAF %d in line 125\n", rank);
        for (int block_idx = 0; block_idx < last_block_idx + di + d + blockSize; block_idx += blockSize)
        {
            // which phase is leaf in?
            else if (block_idx < last_block_idx && block_idx < d + di) // phase 1 - between node's first block arriving at root and receiving first reduced blocks
            {
                MPI_Sendrecv(sendbuf + block_idx, blockSize, MPI_DOUBLE, parent, 0, recvbuf, 0, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else if (block_idx > last_block_idx && block_idx < d + di) // phase 2a - no reduced block incoming yet, but all blocks sent up (for last_block_idx<<d)
            {
                MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, parent, 0, recvbuf, 0, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else if (block_idx >= d + di && block_idx < last_block_idx) // phase 2b - from receiving reduced blocks to sending node's last block upwards;
            {
                MPI_Sendrecv(sendbuf + block_idx, blockSize, MPI_DOUBLE, parent, 0, recvbuf + block_idx - (d + di), blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else if (block_idx >= last_block_idx && block_idx <= d + di + last_block_idx) // phase 3 - starts with all blocks being sent up; ends with last reduced blocks coming in
            {
                MPI_Sendrecv(sendbf, 0, MPI_DOUBLE, parent, 0, recvbuf + block_idx - (d + di), blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        // printf("LEAF %d in line 134\n", rank);
    }

    else
    {
        for (int block_idx = 0; block_idx < last_block_idx + di + d + blockSize; block_idx += blockSize)
        {
            // which phase is middlenode in?
            if (block_idx < d - di || (block_idx > last_block_idx + d - di && block_idx <= d + di)) // phase 1 - until first non-zero block arrives from below || from last upwards send to first reduced receive
            {
                if (childLexists)
                    MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, childL, 0, recvbuf, 0, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (childRexists)
                    MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, childR, 0, recvbuf, 0, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else if (block_idx = d - di) // phase 2: induce offset - receive first block from children
            {
                if (childLexists)
                {
                    MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, childL, 0, tmpbuf, blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + block_idx - (d - di), blockSize, MPI_DOUBLE, MPI_MAX);
                }
                if (childRexists)
                {
                    MPI_Sendrecv(sendbuf, 0, MPI_DOUBLE, childR, 0, tmpbuf, blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + block_idx - (d - di), blockSize, MPI_DOUBLE, MPI_MAX);
                }
            }
            else if (block_idx = d + di + last_block_idx + 1) // phase 8: nothing being received anymore - last send downwards
            {
                if (childLexists)
                {
                    MPI_Sendrecv(recvbuf + last_block_idx, blockSize, MPI_DOUBLE, childL, 0, tmpbuf, 0, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                if (childRexists)
                {
                    MPI_Sendrecv(recvbuf + last_block_idx, blockSize, MPI_DOUBLE, childR, 0, tmpbuf, 0, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            else if (block_idx >= d + di && block_idx <= d + di + last_block_idx && block_idx > d - di + last_block_idx)
            /* Phase 7:
                        1. cond. start receiving stuff from above
                        2. cond. end receiving stuff from above
                        3. cond. nothing being sent up anymore
            */
            {
                if (childLexists)
                {
                    MPI_Sendrecv(recvbuf + block_idx - (d + di), blockSize, MPI_DOUBLE, childL, 0, tmpbuf, 0, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                if (childRexists)
                {
                    MPI_Sendrecv(recvbuf + last_block_idx, blockSize, MPI_DOUBLE, childR, 0, tmpbuf, 0, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

            else if (block_idx > d - di && block_idx <= last_block_idx + (d - di)) // phase 2a - from receiving second block from children to sending last block up
            {
                MPI_Sendrecv(recvbuf + block_idx - (d - di) - blockSize, blockSize, MPI_DOUBLE, parent, 0, tmpbuf, blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (childLexists)
                {
                    MPI_Sendrecv(recvbuf, 0, MPI_DOUBLE, childL, 0, tmpbuf, blockSize, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + block_idx - (d - di), blockSize, MPI_DOUBLE, MPI_MAX);
                }
                if (childRexists)
                {
                    MPI_Sendrecv(recvbuf + block_idx - (d - di) - blockSize, blockSize, MPI_DOUBLE, childR, 0, tmpbuf, blockSize, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + block_idx - (d - di), blockSize, MPI_DOUBLE, MPI_MAX);
                }
            }
            else if (block_idx > d - di && block_idx <= last_block_idx + (d - di)) // phase 2b - receiving first reduced block from parent to receiving last black from children;
            {
            }
            else if (block_idx > d - di && block_idx <= last_block_idx + (d - di)) // phase 2c - receiving zeros from children but still receiving reduced blocks from parent;
            {
            }
            else if (block_idx = last_block_idx + d + blockSize) // redeem offset - send last reduced block to children
            {
                if (childLexists)
                {
                    MPI_Sendrecv(recvbuf + block_idx - d - blockSize, blockSize, MPI_DOUBLE, childL, 0, tmpbuf, 0, MPI_DOUBLE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                if (childRexists)
                {
                    MPI_Sendrecv(recvbuf + block_idx - d - blockSize, blockSize, MPI_DOUBLE, childR, 0, tmpbuf, 0, MPI_DOUBLE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        return MPI_SUCCESS;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int PipelinedTreeduceTreecast(void *sendbuf, void *recvbuf, int count, int blockSize, int size, int rank)
{
    if (size == 1)
    {
        memcpy(recvbuf, sendbuf, count * sizeof(double));
        // could be done more elegantly with pointers.
        // Also the datatype should not be hardcoded
        return MPI_SUCCESS;
    }

    // treeduce(sendbuf, recvbuf, count, blockSize, size, rank);
    return MPI_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    int rank, size, count, blockSize = 10, i;
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
