#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

int main(int argc, char const *argv[])
{
    int dimension, generation, rank, numOfProcesses, blockDimension, i, j, w, z, x, y;
    int right, left, up, down, upRight, upLeft, downRight, downLeft;
    int neighbors;
    int **temp;
    int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3;
    int UP_LEFT = 4, UP_RIGHT = 5, DOWN_LEFT = 6, DOWN_RIGHT = 7;
    double start_time, end_time;
    FILE *inputFile;
    int fileValue;

    MPI_Init(NULL, NULL);
    //get the  number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
    // get the rank 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // every block should have different randomization
    srand(time(NULL) + rank);

    // get the user inputs
    dimension = atoi(argv[1]); // the game of life table is n x n  where n = dimension
    generation = atoi(argv[2]);
    inputFile = fopen(argv[3], "r");

    // each block will be a square table m x m where m = blockDimension
    int sqrtOfProcesses = (int)sqrt(numOfProcesses);
    blockDimension = dimension/sqrtOfProcesses;

    // allocate memory for the current and the new generation block of the current process
    int** currentGen = malloc((blockDimension + 2)*sizeof(int*));
    int** newGen = malloc((blockDimension + 2)*sizeof(int*));
   
    currentGen[0] = malloc((blockDimension + 2)*(blockDimension + 2)*sizeof(int));
    newGen[0] = malloc((blockDimension + 2)*(blockDimension + 2)*sizeof(int));
  
    for (i = 1; i < blockDimension + 2; i++) {
        currentGen[i] = &(currentGen[0][i*(blockDimension + 2)]);
        newGen[i]  = &(newGen[0][i*(blockDimension + 2)]);
    }
    
    // initialize generation from the given input file
    if(inputFile != NULL)
    {
        for (i = 1; i < blockDimension + 1; ++i) {
            for (j = 1; j < blockDimension + 1; ++j) {
                newGen[i][j] = 0;
                 // initialize new generation with zero values
                 fscanf(inputFile, "%d\n", &fileValue);
                 currentGen[i][j] = fileValue;
            }
        }
        fclose(inputFile);
        free(inputFile);
    }
    // initialize generation randomly
    else
    {
        for (i = 1; i < blockDimension + 1; i++) {
            for (j = 1; j < blockDimension + 1; j++) {
                if(rand() % 10 > 5)
                {
                    currentGen[i][j] = 0;
                    newGen[i][j] = 0; 
                }
                else
                {
                    currentGen[i][j] = 1;
                    newGen[i][j] = 1;
                }
            }
        }
    }

    /* START-find the neighbors for the current rank*/
    // right neighbor
    if((rank + 1) % sqrtOfProcesses == 0 )
    {
        right = rank - sqrtOfProcesses + 1;
    }
    else
    {
        right =  rank + 1;
    }
    // left neighbor
    if(rank % sqrtOfProcesses == 0)
    {
        left = rank + sqrtOfProcesses - 1;
    }
    else
    {
        left = rank - 1;
    }
    // up neighbor
    if(rank < sqrtOfProcesses)
    {
        up = rank + numOfProcesses - sqrtOfProcesses;
    }
    else
    {
        up =  rank - sqrtOfProcesses;
    }
    // down neighbor
    if(rank >= numOfProcesses - sqrtOfProcesses)
    {
        down =  rank - numOfProcesses + sqrtOfProcesses;
    }
    else
    {
        down = rank + sqrtOfProcesses;
    }
    // upLeft neighbor
    if (left >= sqrtOfProcesses) 
    {
        upLeft = left - sqrtOfProcesses;
    }
    else 
    {
        upLeft = numOfProcesses - sqrtOfProcesses + left;
    }
    //  upRight neighbor
    if (right >= sqrtOfProcesses)
    {
        upRight = right - sqrtOfProcesses;
    }   
    else 
    {
        upRight = numOfProcesses - sqrtOfProcesses + right;
    }
    // downleft neighbor
    if (left < numOfProcesses - sqrtOfProcesses) 
    {
        downLeft = left + sqrtOfProcesses;
    }
    else 
    {
        downLeft = sqrtOfProcesses - numOfProcesses + left;
    }
    // downRight neighbor
    if (right < numOfProcesses - sqrtOfProcesses) 
    {
        downRight = right + sqrtOfProcesses;
    }
    else 
    {
        downRight = sqrtOfProcesses - numOfProcesses + right;
    }
    /* END-find the neighbors for the current rank*/

    MPI_Datatype columnDatatype; //column data type;
    MPI_Datatype rowDatatype; //row data type;
    
    MPI_Type_vector(blockDimension, 1, blockDimension + 2, MPI_INT, &columnDatatype);
    MPI_Type_commit(&columnDatatype);

    MPI_Type_vector(blockDimension, 1, 1, MPI_INT, &rowDatatype);
    MPI_Type_commit(&rowDatatype);

    /* START - first handle for persistent communication */
    MPI_Request requests_send_firstHandle[8], requests_receive_firstHandle[8];
    MPI_Status status_send_firstHandle[8],status_receive_firstHandle[8];

    // receive data from the neihgbors
    // the tag indicates the receiver's place(down, left, right...) in relation to the sender
    // for example the tag DOWN indicates that the receiver is below the sender
    MPI_Recv_init(&(currentGen[0][1]), 1, rowDatatype, up, DOWN, MPI_COMM_WORLD, &requests_receive_firstHandle[0]);
    MPI_Recv_init(&(currentGen[blockDimension + 1][1]), 1, rowDatatype, down, UP, MPI_COMM_WORLD, &requests_receive_firstHandle[1]);
    MPI_Recv_init(&(currentGen[1][0]), 1, columnDatatype, left, RIGHT, MPI_COMM_WORLD, &requests_receive_firstHandle[2]);
    MPI_Recv_init(&(currentGen[1][blockDimension + 1]), 1, columnDatatype, right, LEFT, MPI_COMM_WORLD, &requests_receive_firstHandle[3]);
    MPI_Recv_init(&(currentGen[0][0]), 1, MPI_INT, upLeft, DOWN_RIGHT, MPI_COMM_WORLD, &requests_receive_firstHandle[4]);
    MPI_Recv_init(&(currentGen[0][blockDimension + 1]), 1, MPI_INT, upRight, DOWN_LEFT, MPI_COMM_WORLD, &requests_receive_firstHandle[5]);
    MPI_Recv_init(&(currentGen[blockDimension + 1][0]), 1, MPI_INT, downLeft, UP_RIGHT, MPI_COMM_WORLD, &requests_receive_firstHandle[6]);
    MPI_Recv_init(&(currentGen[blockDimension + 1][blockDimension + 1]), 1, MPI_INT, downRight, UP_LEFT, MPI_COMM_WORLD, &requests_receive_firstHandle[7]);

    // send data to the neighbors
    MPI_Send_init(&(currentGen[1][1]), 1, rowDatatype, up, UP, MPI_COMM_WORLD, &requests_send_firstHandle[0]);
    MPI_Send_init(&(currentGen[blockDimension][1]), 1, rowDatatype, down, DOWN, MPI_COMM_WORLD, &requests_send_firstHandle[1]);
    MPI_Send_init(&(currentGen[1][1]), 1, columnDatatype, left, LEFT, MPI_COMM_WORLD, &requests_send_firstHandle[2]);
    MPI_Send_init(&(currentGen[1][blockDimension]), 1, columnDatatype, right, RIGHT, MPI_COMM_WORLD, &requests_send_firstHandle[3]);
    MPI_Send_init(&(currentGen[1][1]), 1, MPI_INT, upLeft, UP_LEFT, MPI_COMM_WORLD, &requests_send_firstHandle[4]);
    MPI_Send_init(&(currentGen[blockDimension][1]), 1, MPI_INT, downLeft, DOWN_LEFT, MPI_COMM_WORLD, &requests_send_firstHandle[5]);
    MPI_Send_init(&(currentGen[1][blockDimension]), 1, MPI_INT, upRight, UP_RIGHT, MPI_COMM_WORLD, &requests_send_firstHandle[6]);
    MPI_Send_init(&(currentGen[blockDimension][blockDimension]), 1, MPI_INT, downRight, DOWN_RIGHT, MPI_COMM_WORLD, &requests_send_firstHandle[7]);
    /* END - first handle for persistent communication */

    /* START - second handle for persistent communication */
    MPI_Request requests_send_secondHandle[8], requests_receive_secondHandle[8];
    MPI_Status status_send_secondHandle[8],status_receive_secondHandle[8];

    MPI_Recv_init(&(newGen[0][1]), 1, rowDatatype, up, DOWN, MPI_COMM_WORLD, &requests_receive_secondHandle[0]);
    MPI_Recv_init(&(newGen[blockDimension + 1][1]), 1, rowDatatype, down, UP, MPI_COMM_WORLD, &requests_receive_secondHandle[1]);
    MPI_Recv_init(&(newGen[1][0]), 1, columnDatatype, left, RIGHT, MPI_COMM_WORLD, &requests_receive_secondHandle[2]);
    MPI_Recv_init(&(newGen[1][blockDimension + 1]), 1, columnDatatype, right, LEFT, MPI_COMM_WORLD, &requests_receive_secondHandle[3]);
    MPI_Recv_init(&(newGen[0][0]), 1, MPI_INT, upLeft, DOWN_RIGHT, MPI_COMM_WORLD, &requests_receive_secondHandle[4]);
    MPI_Recv_init(&(newGen[0][blockDimension + 1]), 1, MPI_INT, upRight, DOWN_LEFT, MPI_COMM_WORLD, &requests_receive_secondHandle[5]);
    MPI_Recv_init(&(newGen[blockDimension + 1][0]), 1, MPI_INT, downLeft, UP_RIGHT, MPI_COMM_WORLD, &requests_receive_secondHandle[6]);
    MPI_Recv_init(&(newGen[blockDimension + 1][blockDimension + 1]), 1, MPI_INT, downRight, UP_LEFT, MPI_COMM_WORLD, &requests_receive_secondHandle[7]);

    MPI_Send_init(&(newGen[1][1]), 1, rowDatatype, up, UP, MPI_COMM_WORLD, &requests_send_secondHandle[0]);
    MPI_Send_init(&(newGen[blockDimension][1]), 1, rowDatatype, down, DOWN, MPI_COMM_WORLD, &requests_send_secondHandle[1]);
    MPI_Send_init(&(newGen[1][1]), 1, columnDatatype, left, LEFT, MPI_COMM_WORLD, &requests_send_secondHandle[2]);
    MPI_Send_init(&(newGen[1][blockDimension]), 1, columnDatatype, right, RIGHT, MPI_COMM_WORLD, &requests_send_secondHandle[3]);
    MPI_Send_init(&(newGen[1][1]), 1, MPI_INT, upLeft, UP_LEFT, MPI_COMM_WORLD, &requests_send_secondHandle[4]);
    MPI_Send_init(&(newGen[blockDimension][1]), 1, MPI_INT, downLeft, DOWN_LEFT, MPI_COMM_WORLD, &requests_send_secondHandle[5]);
    MPI_Send_init(&(newGen[1][blockDimension]), 1, MPI_INT, upRight, UP_RIGHT, MPI_COMM_WORLD, &requests_send_secondHandle[6]);
    MPI_Send_init(&(newGen[blockDimension][blockDimension]), 1, MPI_INT, downRight, DOWN_RIGHT, MPI_COMM_WORLD, &requests_send_secondHandle[7]);
    /* END - second handle for persistent communication */
  
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    for(i=0; i < generation; i++)
    {
        if(i % 2 == 0)
        {
            MPI_Startall(8, requests_receive_firstHandle);
            MPI_Startall(8, requests_send_firstHandle);
        }
        else
        {
            MPI_Startall(8, requests_receive_secondHandle);
            MPI_Startall(8, requests_send_secondHandle);
        }

        // calculate the internal cells
        for (x = 2; x < blockDimension; x++) {
            for (y = 2; y < blockDimension; y++) {
                neighbors = 0;
                // count neighbors
                for (w = x-1; w <= x+1; w++) {
                    for (z = y-1; z <= y+1; z++) {
                        if(!(w == x && z == y)){
                            neighbors+=currentGen[w][z];
                        }
                    }
                }

                if (neighbors < 2 || neighbors > 3 || (currentGen[x][y] == 0 && neighbors == 2))
                {
                    newGen[x][y] = 0;
                } 
                else if (neighbors == 3 || (currentGen[x][y] == 1 && neighbors >= 2))
                {
                    newGen[x][y] = 1;
                }
            }
        }

        // wait for all receives to complete
        if(i % 2 == 0)
        {
            MPI_Waitall(8,requests_receive_firstHandle,status_receive_firstHandle);
        }
        else
        {
            MPI_Waitall(8,requests_receive_secondHandle,status_receive_secondHandle);
        }
        
        // calculate the external cells after receiving the required data
        for(x = 1; x < blockDimension + 1; x++){
            // for first row
            neighbors= 0;
            for (w = 0; w <= 2; w++) 
            {
                for (z = x - 1; z <= x + 1; z++) 
                {
                    if(!(w == 1 && z == x))
                    {
                        neighbors += currentGen[w][z];
                    }
                }
            }

            // START - update first row
            if (neighbors < 2 || neighbors > 3 || (currentGen[1][x] == 0 && neighbors == 2)) 
            {
                newGen[1][x] = 0;
            } 
            else if (neighbors == 3 || (currentGen[1][x] == 1 && neighbors >= 2))
            {
                newGen[1][x] = 1;
            }
            // END - update first row

            // for last row
            neighbors = 0;
            for (w = blockDimension-1; w <= blockDimension + 1; w++) 
            {
                for (z = x - 1; z <= x + 1; z++) 
                {
                    if(!(w == blockDimension && z == x))
                    {
                        neighbors += currentGen[w][z];
                    }
                }
            }

            // START - update last row
            if (neighbors < 2 || neighbors > 3 || (currentGen[blockDimension][x] == 0 && neighbors == 2))
            {
                newGen[blockDimension][x] = 0;
            } 
            else if (neighbors == 3 || (currentGen[blockDimension][x] == 1 && neighbors >= 2))
            {
                newGen[blockDimension][x] = 1;
            }
            // END - update last row
        
            // for left column
            neighbors = 0;
            for (w = x-1; w <= x+1; w++) 
            {
                for (z = 0; z <= 2; z++) 
                {
                    if(!(w == x && z == 1))
                    {
                        neighbors  += currentGen[w][z];
                    }
                }
            }

            // START - update left column
            if (neighbors < 2 || neighbors > 3 || (currentGen[x][1] == 0 && neighbors == 2))
            {
                newGen[x][1] = 0;
            } 
            else if (neighbors == 3 || (currentGen[x][1] == 1 && neighbors >= 2))
            {
                newGen[x][1] = 1;
            }
            // END - update left column

            // for right column
            neighbors = 0;
            for (w = x-1; w <= x+1; w++) 
            {
                for (z = blockDimension - 1; z <= blockDimension + 1; z++) 
                {
                    if(!(w == x && z == blockDimension))
                    {
                        neighbors += currentGen[w][z];
                    }
                }
            }

            // START - update right column
            if (neighbors < 2 || neighbors > 3 || (currentGen[x][blockDimension] == 0 && neighbors == 2))
            {
                newGen[x][blockDimension] = 0;
            } 
            else if (neighbors == 3 || (currentGen[x][blockDimension] == 1 && neighbors >= 2))
            {
                newGen[x][blockDimension] = 1;
            }
            // END - update right column
        }

        // make the new generation the current generation
        temp = currentGen;
        currentGen = newGen;
        newGen = temp;

        if(i % 2 == 0)
        {
            MPI_Waitall(8, requests_send_firstHandle, status_send_firstHandle);
        }
        else
        {
            MPI_Waitall(8, requests_send_secondHandle, status_send_secondHandle);

        }

    }

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    
    if (rank == 0) {
        printf("Execution Time: %f\n", end_time - start_time);
    }

    // free current and new generation table
    free(currentGen[0]);
    free(newGen[0]);
    free(currentGen);
    free(newGen);

    MPI_Type_free(&columnDatatype);
    MPI_Type_free(&rowDatatype);

    for(int r = 0; r < 8 ; r++)
    {
        MPI_Request_free (&requests_receive_firstHandle[r]);
        MPI_Request_free (&requests_receive_secondHandle[r]);
        MPI_Request_free (&requests_send_firstHandle[r]);
        MPI_Request_free (&requests_send_secondHandle[r]);
    }

    MPI_Finalize();
    return 0;
}