///////////////
// DATA TYPE //
///////////////

#define TUW_TYPE MPI_FLOAT
typedef float tuwtype_t;

////////////////
// Exercise 1 //
////////////////

int allReduce1(tuwtype_t *sendbuf, int sendcount, tuwtype_t *recvbuf, int recvcount, MPI_Comm comm)
{
    MPI_Reduce(sendbuf, recvbuf, sendcount, TUW_TYPE, MPI_MAX, 0, comm);
    MPI_Bcast(recvbuf, sendcount, TUW_TYPE, 0, comm);
    return MPI_SUCCESS;
}

////////////////
// Exercise 2 //
////////////////

// reduce method
int reduce(tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int trunc_chunk_idx = floor(count / blockSize) * blockSize;
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbuf_trunc = (tuwtype_t *)malloc(count % blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;

    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(tmpbuf, blockSize, TUW_TYPE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, TUW_TYPE, MPI_MAX);
        }
        if (count % blockSize != 0)
        {
            MPI_Recv(tmpbuf_trunc, count % blockSize, TUW_TYPE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, MPI_MAX);
        }
    }
    else if (rank == size - 1)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(recvbuf + b, blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (count % blockSize != 0)
        {
            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        if (count >= blockSize)
        {
            MPI_Recv(tmpbuf, blockSize, TUW_TYPE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Reduce_local(tmpbuf, recvbuf, blockSize, TUW_TYPE, MPI_MAX);

            for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
            {
                MPI_Sendrecv(recvbuf + b - blockSize, blockSize, TUW_TYPE, rank - 1, 0, tmpbuf, blockSize, TUW_TYPE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, TUW_TYPE, MPI_MAX);
            }
            MPI_Send(recvbuf + trunc_chunk_idx - blockSize, blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (count % blockSize != 0)
        {
            MPI_Recv(tmpbuf_trunc, count % blockSize, TUW_TYPE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, MPI_MAX);

            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD);
        }
    }
    // free temporary memory
    free(tmpbuf);
    free(tmpbuf_trunc);
    return MPI_SUCCESS;
}

// broadcast method
int bcast(tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int trunc_chunk_idx = floor(count / blockSize) * blockSize;
    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(recvbuf + b, blockSize, TUW_TYPE, rank + 1, 0, MPI_COMM_WORLD);
        }
        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, rank + 1, 0, MPI_COMM_WORLD);
        }
    }
    else if (rank == size - 1)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(recvbuf + b, blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Recv(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        if (count >= blockSize)
        {
            MPI_Recv(recvbuf, blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int b = blockSize; b < trunc_chunk_idx; b += blockSize)
            {
                MPI_Sendrecv(recvbuf + b - blockSize, blockSize, TUW_TYPE, rank + 1, 0, recvbuf + b, blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            // send last 'complete block'
            MPI_Send(recvbuf + trunc_chunk_idx - blockSize, blockSize, TUW_TYPE, rank + 1, 0, MPI_COMM_WORLD);
        }
        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Recv(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    return MPI_SUCCESS;
}

int allReduce2(tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{

    // account for trivial case
    if (size == 1)
    {
        return MPI_SUCCESS;
    }

    reduce(recvbuf, count, blockSize, size, rank);
    bcast(recvbuf, count, blockSize, size, rank);

    return MPI_SUCCESS;
}

////////////////
// exercise 3 //
////////////////

int allReduce3(tuwtype_t *sendbuf, int count, int blockSize, int size, int rank)
{
    if (size == 1)
    {
        return MPI_SUCCESS;
    }

    int b = ceil(((double)count) / blockSize);
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));

    int child = rank + 1;
    int parent = rank - 1;
    int d_max = size;
    int di = rank;

    int childExists = (child < size);

    // handle receiving of zeros
    MPI_Status status;
    int elements;

    for (int j = 0; j <= b + d_max; ++j) // NOTE: j runs to b + d_max to match sends and receives
    {
        if (di < d_max) // non-leaf nodes
        {
            if (((j - (di + 1)) < 0) || ((j - (di + 1)) >= b)) // 'out-of-bounds' condition (only affects sends)
            {
                if (childExists)
                {
                    MPI_Sendrecv(sendbuf, 0, TUW_TYPE, child, 0, tmpbuf, blockSize, TUW_TYPE, child, 0, MPI_COMM_WORLD, &status);
                    if (j < b)
                    {
                        MPI_Get_count(&status, TUW_TYPE, &elements);
                        MPI_Reduce_local(tmpbuf, sendbuf + j * blockSize, elements, TUW_TYPE, MPI_MAX);
                    }
                }
            }
            else
            {
                if (childExists)
                {
                    MPI_Sendrecv(sendbuf + (j - (di + 1)) * blockSize, blockSize, TUW_TYPE, child, 0, tmpbuf, blockSize, TUW_TYPE, child, 0, MPI_COMM_WORLD, &status);
                    if (j < b)
                    {
                        MPI_Get_count(&status, TUW_TYPE, &elements);
                        MPI_Reduce_local(tmpbuf, sendbuf + j * blockSize, elements, TUW_TYPE, MPI_MAX);
                    }
                }
            }
        }
        if (di > 0) // non-root nodes
        {
            if (j > b + di) // 'do-nothing' condition since j runs to d_max in our code
            {
                MPI_Sendrecv(sendbuf, 0, TUW_TYPE, parent, 0, tmpbuf, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, &status);
            }
            else
            {
                if (j >= b)
                {
                    MPI_Sendrecv(sendbuf, 0, TUW_TYPE, parent, 0, tmpbuf, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, &status);
                }
                else
                {
                    MPI_Sendrecv(sendbuf + j * blockSize, blockSize, TUW_TYPE, parent, 0, tmpbuf, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, &status);
                }
                if ((j - di) >= 0 && (j - di) < b)
                {
                    MPI_Get_count(&status, TUW_TYPE, &elements);
                    MPI_Reduce_local(tmpbuf, sendbuf + (j - di) * blockSize, elements, TUW_TYPE, MPI_MAX);
                }
            }
        }
    }
    free(tmpbuf);
    return MPI_SUCCESS;
}

////////////////
// exercise 4 //
////////////////

int treeduce(tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int childL = 2 * rank + 1;
    int childR = 2 * rank + 2;
    int childLexists = (childL < size);
    int childRexists = (childR < size);

    tuwtype_t d = ceil(log2(size + 1) - 1); // binary tree depth
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

    int trunc_chunk_idx = floor(count / blockSize) * blockSize; // index of first truncated chunk
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));
    tuwtype_t *tmpbuf_trunc = (tuwtype_t *)malloc(count % blockSize * sizeof(tuwtype_t));
    tuwtype_t tmp;
    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            if (childLexists)
            {
                MPI_Recv(tmpbuf, blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, TUW_TYPE, MPI_MAX);
            }
            if (childRexists)
            {
                MPI_Recv(tmpbuf, blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, TUW_TYPE, MPI_MAX);
            }
        }
        if (count % blockSize != 0)
        {
            if (childLexists)
            {
                MPI_Recv(tmpbuf_trunc, count % blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, MPI_MAX);
            }
            if (childRexists)
            {
                MPI_Recv(tmpbuf_trunc, count % blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, MPI_MAX);
            }
        }
    }
    else if ((!childRexists) && (!childLexists)) // check if this is a leaf-node
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Send(recvbuf + b, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD);
        }
        if (count % blockSize != 0)
        {
            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        if (count >= blockSize)
        {
            for (int b = 0; b < trunc_chunk_idx; b += blockSize)
            {
                if (childLexists)
                {
                    MPI_Recv(tmpbuf, blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, TUW_TYPE, MPI_MAX);
                }
                if (childRexists)
                {
                    MPI_Recv(tmpbuf, blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Reduce_local(tmpbuf, recvbuf + b, blockSize, TUW_TYPE, MPI_MAX);
                }
                MPI_Send(recvbuf + b, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD);
            }
        }
        if (count % blockSize != 0)
        {
            if (childLexists)
            {
                MPI_Recv(tmpbuf_trunc, count % blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, MPI_MAX);
            }
            if (childRexists)
            {
                MPI_Recv(tmpbuf_trunc, count % blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Reduce_local(tmpbuf_trunc, recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, MPI_MAX);
            }
            MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD);
        }
    }
    free(tmpbuf);
    free(tmpbuf_trunc);
    return MPI_SUCCESS;
}

// binary tree broadcast
int treecast(tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    int childL = 2 * rank + 1;
    int childR = 2 * rank + 2;
    int childLexists = (childL < size);
    int childRexists = (childR < size);

    tuwtype_t d = ceil(log2(size + 1) - 1); // binary tree depth
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

    int trunc_chunk_idx = floor(count / blockSize) * blockSize;

    if (rank == 0)
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            if (childLexists)
                MPI_Send(recvbuf + b, blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD);
            if (childRexists)
                MPI_Send(recvbuf + b, blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD);
        }

        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            if (childLexists)
                MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD);
            if (childRexists)
                MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD);
        }
    }

    else if ((!childRexists) && (!childLexists)) // check if this is a leaf-node
    {
        for (int b = 0; b < trunc_chunk_idx; b += blockSize)
        {
            MPI_Recv(recvbuf + b, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Recv(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        if (count >= blockSize)
        {
            for (int b = 0; b < trunc_chunk_idx; b += blockSize)
            {
                MPI_Recv(recvbuf + b, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (childLexists)
                    MPI_Send(recvbuf + b, blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD);
                if (childRexists)
                {
                    MPI_Send(recvbuf + b, blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD);
                }
            }
        }
        /* added this to account for blockSizes that don't divide count */
        if (count % blockSize != 0)
        {
            MPI_Recv(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (childLexists)
                MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD);
            if (childRexists)
                MPI_Send(recvbuf + trunc_chunk_idx, count % blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD);
        }
    }
    return MPI_SUCCESS;
}

int allReduce4(tuwtype_t *recvbuf, int count, int blockSize, int size, int rank)
{
    if (size == 1)
    {
        return MPI_SUCCESS;
    }

    treeduce(recvbuf, count, blockSize, size, rank);
    treecast(recvbuf, count, blockSize, size, rank);

    return MPI_SUCCESS;
}

////////////////
// exercise 5 //
////////////////

int allReduce5(tuwtype_t *sendbuf, int count, int blockSize, int size, int rank)
{
    if (size == 1)
    {
        return MPI_SUCCESS;
    }

    int b = ceil(((double)count) / blockSize);
    tuwtype_t *tmpbuf = (tuwtype_t *)malloc(blockSize * sizeof(tuwtype_t));

    int childL = 2 * rank + 1;
    int childR = 2 * rank + 2;

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

    int d_max = ceil(log2(size + 1) - 1); // binary tree depth
    int di;

    // find node depth
    tuwtype_t size_run = 0;
    for (int d_run = 0; d_run < d_max + 1; d_run++)
    {
        size_run += pow(2, d_run);
        if (rank < size_run)
        {
            di = d_run;
            break;
        }
    }

    int childLexists = (childL < size);
    int childRexists = (childR < size);

    // handle receiving of zeros
    MPI_Status status;
    int elements;

    for (int j = 0; j <= b + d_max; ++j) // NOTE: j runs to b + d_max to match sends and receives
    {
        if (di < d_max)
        { // non-leaf nodes
            if (((j - (di + 1)) < 0) || ((j - (di + 1)) >= b))
            { // 'out-of-bounds' condition (affects send only)
                if (childLexists)
                {
                    MPI_Sendrecv(sendbuf, 0, TUW_TYPE, childL, 0, tmpbuf, blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD, &status);
                    if (j < b)
                    {
                        MPI_Get_count(&status, TUW_TYPE, &elements);
                        MPI_Reduce_local(tmpbuf, sendbuf + j * blockSize, elements, TUW_TYPE, MPI_MAX);
                    }
                }
                if (childRexists)
                {
                    MPI_Sendrecv(sendbuf, 0, TUW_TYPE, childR, 0, tmpbuf, blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD, &status);
                    if (j < b)
                    {
                        MPI_Get_count(&status, TUW_TYPE, &elements);
                        MPI_Reduce_local(tmpbuf, sendbuf + j * blockSize, elements, TUW_TYPE, MPI_MAX);
                    }
                }
            }
            else
            {
                if (childLexists)
                {
                    MPI_Sendrecv(sendbuf + (j - (di + 1)) * blockSize, blockSize, TUW_TYPE, childL, 0, tmpbuf, blockSize, TUW_TYPE, childL, 0, MPI_COMM_WORLD, &status);
                    if (j < b)
                    {
                        MPI_Get_count(&status, TUW_TYPE, &elements);
                        MPI_Reduce_local(tmpbuf, sendbuf + j * blockSize, elements, TUW_TYPE, MPI_MAX);
                    }
                }
                if (childRexists)
                {
                    MPI_Sendrecv(sendbuf + (j - (di + 1)) * blockSize, blockSize, TUW_TYPE, childR, 0, tmpbuf, blockSize, TUW_TYPE, childR, 0, MPI_COMM_WORLD, &status);
                    if (j < b)
                    {
                        MPI_Get_count(&status, TUW_TYPE, &elements);
                        MPI_Reduce_local(tmpbuf, sendbuf + j * blockSize, elements, TUW_TYPE, MPI_MAX);
                    }
                }
            }
        }
        if (di > 0)
        { // non-root nodes
            if (j > b + di)
            { // 'do-nothing' condition since j runs to d_max in our code
                MPI_Sendrecv(sendbuf, 0, TUW_TYPE, parent, 0, tmpbuf, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, &status);
            }
            else
            {
                if (j >= b)
                {
                    MPI_Sendrecv(sendbuf, 0, TUW_TYPE, parent, 0, tmpbuf, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, &status);
                }
                else
                {
                    MPI_Sendrecv(sendbuf + j * blockSize, blockSize, TUW_TYPE, parent, 0, tmpbuf, blockSize, TUW_TYPE, parent, 0, MPI_COMM_WORLD, &status);
                }
                if ((j - di) >= 0 && (j - di) < b)
                {
                    MPI_Get_count(&status, TUW_TYPE, &elements);
                    MPI_Reduce_local(tmpbuf, sendbuf + (j - di) * blockSize, elements, TUW_TYPE, MPI_MAX);
                }
            }
        }
    }
    free(tmpbuf);

    return MPI_SUCCESS;
}