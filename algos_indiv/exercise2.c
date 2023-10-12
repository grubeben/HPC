#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include <assert.h>

#include <mpi.h>

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

// reduce method
int reduce(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int trunc_chunk_idx = floor(count / blockSize) * blockSize; // index of first truncated chunk.. yes funny.
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbuf_trunc = (tuwtype_t *)malloc(count % blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;

    if (rank == 0)
    {
        // printf("start root\n");
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {

            MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            /*
            for (int j = 0; j < blockSize; j++)
            {
                recvbuf[b + j] = tmpbuf[j] > recvbuf[b + j] ? tmpbuf[j] : recvbuf[b + j];
            }
            */
            MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, MPI_DOUBLE, MPI_MAX);
            // printf("root even part done\n");
        }
        if (count % blockSize != 0)
        {

            MPI_Recv(tmpbuf_trunc, count % blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            /*
            for (int j = 0; j < count % blockSize; j++)
            {
                recvbuf[trunc_chunk_idx + j] = tmpbuf_trunc[j] > recvbuf[trunc_chunk_idx + j] ? tmpbuf_trunc[j] : recvbuf[trunc_chunk_idx + j];
            }
            */
            MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, MPI_MAX);
        }
        // printf("root done\n");
    }
    else if (rank == size - 1)
    {
        // printf("send only node start\n");
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(sendbuf + b, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (count % blockSize != 0)
        {
            MPI_Send(sendbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        // printf("send only node done\n");
    }
    else
    {
        if (count >= blockSize)
        {
            // printf("regular node %f start\n", rank);
            MPI_Recv(tmpbuf, blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            /*
            for (int j = 0; j < blockSize; j++)
            {
                recvbuf[j] = tmpbuf[j] > recvbuf[j] ? tmpbuf[j] : recvbuf[j];
            }
            */
            MPI_Reduce_local(tmpbuf, recvbuf, blockSize, MPI_DOUBLE, MPI_MAX);

            for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
            {

                MPI_Sendrecv(recvbuf + b - blockSize, blockSize, MPI_DOUBLE, rank - 1, 0, tmpbuf, blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                /*
                for (int j = 0; j < blockSize; j++)
                {
                    recvbuf[b + j] = tmpbuf[j] > recvbuf[b + j] ? tmpbuf[j] : recvbuf[b + j];
                }
                */
                MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, MPI_DOUBLE, MPI_MAX);
            }
            // printf("regular node %f even part done\n", rank);
            MPI_Send(recvbuf + trunc_chunk_idx - blockSize, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (count % blockSize != 0)
        {
            MPI_Recv(tmpbuf_trunc, count % blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            /*
            for (int j = 0; j < count % blockSize; j++)
            {
                recvbuf[trunc_chunk_idx + j] = tmpbuf_trunc[j] > recvbuf[trunc_chunk_idx + j] ? tmpbuf_trunc[j] : recvbuf[trunc_chunk_idx + j];
            }
            */
            MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, MPI_MAX);

            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            // printf("regular node %f done\n", rank);
        }
    }
    // free temporary memory
    // free(tmpbuf);
    // free(tmpbuf_trunc);
    return MPI_SUCCESS;
}

// broadcast
int bcast(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int trunc_chunk_idx = floor(count / blockSize) * blockSize;

    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(recvbuf + b, blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }

        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
    }
    else if (rank == size - 1)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(recvbuf + b, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Recv(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        if (count >= blockSize)
        {
            MPI_Recv(recvbuf, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("%d\n",rank);
            for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
            {
                // printf("before sendrc%d\n",rank);
                MPI_Sendrecv(recvbuf + b - blockSize, blockSize, MPI_DOUBLE, rank + 1, 0, recvbuf + b, blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // printf("after sendrc%d\n",rank);
            }
            // send last 'complete block'
            MPI_Send(recvbuf + trunc_chunk_idx - blockSize, blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Recv(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    return MPI_SUCCESS;
}

int PipelinedReduceBcast(void *sendbuf, void *recvbuf, int count, int blockSize, int size, int rank)
{
    if (size == 1)
    {
         memcpy(recvbuf, sendbuf, count * sizeof(double)); 
         // could be done more elegantly with pointers. 
         // Also the datatype should not be hardcoded
        return MPI_SUCCESS;
    }

    reduce(sendbuf, recvbuf, count, blockSize, size, rank);
    bcast(sendbuf, recvbuf, count, blockSize, size, rank);

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
    PipelinedReduceBcast(sendbuf, testbuf, count, blockSize, size, rank);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //testbuf = sendbuf;

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

    // free(sendbuf);
    // free(recvbuf);
    // free(testbuf);

    return 0;
}
