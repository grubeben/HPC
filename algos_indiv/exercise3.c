
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include <assert.h>

#include <mpi.h>

#define TUW_TYPE MPI_DOUBLE
typedef double tuwtype_t;

int PipelinedTree_withoutbranches_duceTrunkcast(tuwtype_t *sendbuf, tuwtype_t *recvbuf, int count, int blockSize, int size, int rank) {

    int b = count/blockSize;
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *zerobuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));

    for (int i = 0; i < blockSize; i++) {
        zerobuf[i] = (tuwtype_t)0.0;
    }

    int child = rank + 1;
    int parent=rank-1;

    int d_max = size;
    int di=rank;

    int childExists = (child < size);

    for (int j = 0; j <= b + d_max; ++j) // NOTE: j runs to b + d_max to match sends and receives
    {
        if (di < d_max) { // non-leaf nodes
            if ( ((j - (di + 1)) < 0) || ((j - (di + 1)) >= b) ) { // 'out-of-bounds' BETRIFFT NUR SEND
                if (childExists) {
                    //printf("j=%d: node %d to/from node %d in line 65\n", j, rank, child);
                    MPI_Sendrecv(zerobuf, blockSize, MPI_DOUBLE, child, 0, tmpbuf, blockSize, MPI_DOUBLE, child, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (j < b)
                        MPI_Reduce_local(tmpbuf, sendbuf + j*blockSize, blockSize, MPI_DOUBLE, MPI_MAX);
                }
            } else {
                if (childExists) {
                    //printf("j=%d: node %d to/from node %d in line 78\n", j, rank, child);
                    MPI_Sendrecv(sendbuf+(j-(di+1))*blockSize, blockSize, MPI_DOUBLE, child, 0, tmpbuf, blockSize, MPI_DOUBLE, child, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (j < b)
                        MPI_Reduce_local(tmpbuf, sendbuf + j*blockSize, blockSize, MPI_DOUBLE, MPI_MAX);
                }
            }
        }
        if (di > 0) { // non-root nodes
            if (j > b + di) { // 'do-nothing' condition since j runs to d_max in our code

                    //printf("j=%d: node %d to/from node %d in line 94\n", j, rank, child);
                    MPI_Sendrecv(zerobuf, blockSize, MPI_DOUBLE, parent, 0, tmpbuf, blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            } else {
                //printf("j=%d: node %d sends senfbuf[%d] to node %d\n", j, rank, j*blockSize, parent);
                if (j >= b) {
                    //printf("j=%d: node %d to/from node %d in line 100\n", j, rank, parent);
                    MPI_Sendrecv(zerobuf, blockSize, MPI_DOUBLE, parent, 0, tmpbuf, blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    //printf("j=%d: node %d to/from node %d in line 103\n", j, rank, parent);
                    MPI_Sendrecv(sendbuf + j*blockSize, blockSize, MPI_DOUBLE, parent, 0, tmpbuf, blockSize, MPI_DOUBLE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                if ((j - di) >= 0 && (j - di) < b)
                    //printf("j=%d: node %d to/from node %d in line 103\n", j, rank, parent);
                    MPI_Reduce_local(tmpbuf, sendbuf+(j-di)*blockSize, blockSize, MPI_DOUBLE, MPI_MAX);
            }
        }
    }



/*
    for (int j = 0; j < b + di; ++j)
    {
        if (not leaf) {
            if ((j - (di + 1) < 0) && (j - (di + 1)) > b) {
                MPI_Sendrecv(0, child, t, child)
                MPI_Reduce_local(Y[j], Y[j], t)
                MPI_Sendrecv(0, childR, t, childR)
                MPI_Reduce_local(Y[j], Y[j], t)
            } else {
                MPI_Sendrecv(Y[j-(di+1)], child, t child)
                MPI_Reduce_local(Y[j], Y[j], t)
                MPI_Sendrecv(Y[j-(di+1)], childR, t, childR)
                MPI_Reduce_local(Y[j], Y[j], t)
            }
        }

        if (not root) {
            if ((j - di) < 0 && (j - di) > b) {
                MPI_Sendrecv(Y[j], parent, 0, parent)
            } else {
                MPI_Sendrecv(Y[j], parent, Y[j-di], parent)
            }
        }
    }
    */
}

int main(int argc, char *argv[])
{
    int rank, size, count, blockSize = 4, i;
    tuwtype_t *sendbuf, *recvbuf, *testbuf, *bothbuf;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of elements in buffers
    count = 16;

    sendbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(sendbuf != NULL);
    recvbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(recvbuf != NULL);
    testbuf = (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(testbuf != NULL);
    bothbuf= (tuwtype_t *)malloc(count * sizeof(tuwtype_t));
    assert(bothbuf != NULL);

    srand(time(NULL));



    for (i = 0; i < count; i++) {
        sendbuf[i] = (tuwtype_t)rank;
        bothbuf[i] = (tuwtype_t)rank;
    }
    for (i = 0; i < count; i++)
        recvbuf[i] = (tuwtype_t)-1;
    for (i = 0; i < count; i++)
        testbuf[i] = (tuwtype_t)-1;

    PipelinedTree_withoutbranches_duceTrunkcast(bothbuf, testbuf, count, blockSize, size, rank);
    MPI_Allreduce(sendbuf, recvbuf, count, TUW_TYPE, MPI_MAX, MPI_COMM_WORLD);



    printf("reference  custom\n");
    for (i = 0; i < count; i++)
    {
        // assert(recvbuf[i] == testbuf[i]);
        printf("%f, %f\n", recvbuf[i], bothbuf[i]);
    }
    printf("\n");



    

    MPI_Finalize();

    free(bothbuf);
    free(sendbuf);
    free(recvbuf);
    free(testbuf);

    return 0;
}