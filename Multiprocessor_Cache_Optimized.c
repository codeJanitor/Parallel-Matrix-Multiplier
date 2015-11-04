/**
 * @author Caleb Mays and Nathan Brown
 * @date 16 Nov 2010
 * Description: Parallel - project 4, part 3. Matrix Multiplication.
 *              Provides a multiprocessor implementation of matrix multiplication,
 *              which recursively splits each sub-matrix 
 *              to optimize local caching performance.
 *              Since this code focuses on performance, some sacrifices have
 *              been made in terms of readability and modularity.
 *              Add one more matrix allocation, 
 *              and subroutines will certainly be required!
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <sys/stat.h>
#include <mpi.h>


#define DATA_MSG 0
#define PROMPT_MSG 1
#define RESPONSE_MSG 2
#define BUFFER_SIZE 10

#define BLOCK_LOW(id, p, n) ((id) * (n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)


void MatrixMultiply(int l, int m, int n,
                    int **cMatrix, int **aMatrix, int **bMatrix,
                    int crow, int ccol,
                    int arow, int acol,
                    int brow, int bcol,
                    int id);

void print_full_matrix(int **array, int row, MPI_Comm comm);
void print_matrix(int **array, int low_bound, int rows, int cols);


int main (int argc, char *argv[]) {
    
	// Columns of A must equal rows of B.
	// C will have rowsA x colsB.
	int *rowsA, *colsA, *rowsB, *colsB, *rowsC, *colsC, *n;
	int     **a;     
	int     *aStorage;
	int 	**b;
	int		*bStorage;
	int 	**c;
	int 	*cStorage;
    int     **temp;
    int     *tempStorage;
    int     local_rows, local_cols;
	int i, j, k;
    
    MPI_Status status, status2;
    MPI_Request handle, handle2;
    
	if (argc != 2){
             printf("not correct set of arguments");
             return 0;
    }
	char * inputFileName =  argv[1];
	FILE * inputFile;
	int id;
	int p;
	double elapsed_time;
	
	
	MPI_Init(&argc, &argv);	
	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
    char numBuffer[BUFFER_SIZE];
    
	// Get the memory I need.
    n = malloc(sizeof(int));
	rowsA = (int*)malloc(sizeof(int));
	colsA = (int*)malloc(sizeof(int));
	rowsB = (int*)malloc(sizeof(int));
	colsB = (int*)malloc(sizeof(int));	
	rowsC = (int*)malloc(sizeof(int));
	colsC = (int*)malloc(sizeof(int));
	
	if (id == (p-1)) {
		inputFile = fopen (inputFileName, "r");
        
		if (inputFile == NULL) {
				   printf("No file pointer!\n");
				   return 0;
		} else {
			// This is going to be all square. also, a power of two.
            // rowsC will be the size for everybody else to know about.
            *n = atoi(fgets(numBufffer, BUFFER_SIZE, inputFile));
		}
	}
    
    MPI_Bcast(n, 1, MPI_INT, p - 1, MPI_COMM_WORLD);
    
    local_rows = BLOCK_SIZE(id, p, *n);
    local_cols = *n;
    
    // So now everybody knows exactly what size each matrix is
    // that they should be dealing with.
    // Everyone should have the same stuff for A, B, and C matrices.
    *colsB = *colsC = *colsA = local_cols;
    *rowsB = *rowsA = *rowsC = local_rows;

	// All of the memory allocation business and addressing business.
	a = (int**)malloc((*rowsA) * sizeof(int*));
	b = (int**)malloc((*rowsB) * sizeof(int*));
	c = (int**)malloc((*rowsC) * sizeof(int*));
    temp = (int**)malloc((*rowsC) * sizeof(int*));
    
	aStorage = (int*)malloc((*rowsA) * (*colsA) * sizeof(int));
	bStorage = (int*)malloc((*rowsB) * (*colsB) * sizeof(int));
	cStorage = (int*)malloc((*rowsC) * (*colsC) * sizeof(int));
    tempStorage = (int*)malloc((*rowsC) * (*colsC) * sizeof(int));
    
	for (i = 0; i < (*rowsA); i++) {
		a[i] = &aStorage[i * (*colsA)];
	}
    
	for (i = 0; i < (*rowsB); i++) {
		b[i] = &bStorage[i * (*colsB)];
        temp[i] = &tempStorage[i * (*colsB)];
	}
    
	for (i = 0; i < (*rowsC); i++) {
		c[i] = &cStorage[i * (*colsC)];
	}	
	

	// Reading in and assigning the array.
    if (id = p - 1) {
        for (i = 0; i < p - 1; i++) {
            for(j = 0; j < local_rows; j++) {
                for(k = 0; k < local_cols; k++) {
                    a[j][k] = fgetc(inputFile) - 48;
                }
                fgetc(inputFile);
                // Assume only single-digit values appear in the matrix file.
            }
            // We have a block of data, and I want to send it out
            // to whomever cares about it.
            MPI_Send(*a, local_rows * local_cols,
                     MPI_INT, i, DATA_MSG, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(*a, local_rows * local_cols, MPI_INT, p - 1, DATA_MSG, MPI_COMM_WORLD, &status);
    }
    
    // We're going to multiply a matrix against itself, so B will equal A.
    // Also, zero out C.
    for(i = 0; i < local_rows, i++) {
        for (j = 0; j < local_cols; j++) {
            b[i][j] = a[i][j];
            c[i][j] = 0;
        }
    }

	// Time it, do it, stop timing.
    MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = -MPI_Wtime();
    
    // According to the lecture, iterate (p - 1) times.
    for (i = 0; i < p - 1; i++) {
        // MPI_Isend b to P + 1 (P0 if last)...a mod of some sort??
        // calculate partial C on the old stuff.
        MPI_Isend(*b, local_rows * local_cols,
                  MPI_INT, (id + 1) % p, DATA_MSG, MPI_COMM_WORLD, &handle2);
        
        // MPI_Irecv b from P - 1 (or last if P0)
        // into a temporary overlap buffer.
        MPI_Irecv(*temp, local_rows * local_cols,
                  MPI_INT, (id + p - 1) % p, DATA_MSG, MPI_COMM_WORLD, &handle);
        
        // Calculate partial C on the old stuff.
        MatrixMultiply(*rowsA, *rowsB, *colsB,
                       c, a, b,
                       0, 0, 0,
                       (*n/p) * ((id - i + p) % p),
                       0, 0, id);
        
        // Wait.
        MPI_Wait(&handle2, &status2);
        MPI_Wait(&handle, &status);
        
        for (k = 0; k < local_rows; k++) {
            for(j = 0; j < local_cols; j++) {
                b[k][j] = temp[k][j];
            }
        }
    }
    
    // Calculate the last partial C on the final stuff.
    MatrixMultiply(*rowsA, *rowsB, *colsB,
                   c, a, b,
                   0, 0, 0,
                   (*n/p) * ((id - i + p) % p),
                   0, 0, id);
    
    elapsed_time += MPI_Wtime();
    
    // Get the information we need.
    
	if (id == 0) {
		FILE *file, *file2;
		file = fopen("myAnswers.csv", "a+");
        file2 = fopen("myAnswers2.csv", "a+");
		fprintf(file, "%d, %d, %10.6f\n", p, *n, elapsed_time);
   		fprintf(file2, "%d, %d, %10.6f\n", *n, p, elapsed_time);
		fclose(file);
        fclose(file2);
	}
	
    // Free memory and exit.
	free(a);
	free(aStorage);
	free(b);
	free(bStorage);
	free(c);
	free(cStorage);
    free(temp);
    free(tempStorage);
	free(rowsA);
	free(rowsB);
	free(rowsC);
	free(colsA);
	free(colsB);
	free(colsC);
	MPI_Finalize();
    
	return 0;
}

/*
void MatrixMultiply(int l, int m, int n, int **cMatrix, int **aMatrix, int **bMatrix, int Rowc,
     int Colc, int Rowa, int Cola, int Rowb, int Colb){
    */
void MatrixMultiply(int l, int m, int n,
                    int **cMatrix, int **aMatrix, int **bMatrix,
                    int crow, int ccol,
                    int arow, int acol,
                    int brow, int bcol,
                    int id) {
    
	int i, j, k;
    int lhalf[3];
    int mhalf[3];
    int nhalf[3];
    
    lhalf[0] = 0; lhalf[1] = l/2; lhalf[2] = l - l/2;
    mhalf[0] = 0; mhalf[1] = m/2; mhalf[2] = m - m/2;
    nhalf[0] = 0; nhalf[1] = n/2; nhalf[2] = n - n/2;
    
    // If my B portion is larger than 256,
    // split it up again with a recursive call.
    if (n * m > 256) {
        for (i = 0; i < 2; i++) {
            for(j = 0; j < 2; j++) {
                for(k = 0; k < 2; k++) {
                    MatrixMultiply(lhalf[i + 1], mhalf[k + 1], nhalf[j + 1],
                                   cMatrix, aMatrix, bMatrix,
                                   crow + lhalf[i], ccol + nhalf[j],
                                   arow + lhalf[j], acol+mhalf[k],
                                   brow + mhalf[k], bcol + nhalf[j],
                                   id);
                }
            }
        }
        
    } else {
        for(i = 0; i < l; i++) {
            for(j = 0; j < n; j++) {
                for(k = 0; k < m; k++) {
                    cMatrix[i + crow][j + ccol] += aMatrix[i + arow][k + acol] * bMatrix[k + brow][j + bcol];
                }
            }
        }
    }
    
    return;
}


/*
print_full_matrix
 parameters:
    an array.
    a number of rows (for the size. We get the number of columns later.
    MPI_COMM_WORLD
 */
void print_full_matrix(int** array, int rows, MPI_Comm comm) {
    int** buffer;
    int* buffStorage;
    int i, j, k, x;
    int prompt;
    int localrows, local_cols;
    int p, id;
    MPI_Status status;
    
    MPI_Comm_rank(comm, &id);
    MPI_Comm_size(comm, &p);
    
    local_rows = BLOCK_SIZE(id, p, rows);
    local_cols = rows;
    
    if (id == 0) {
        printf("Processor %d: Printing Final Matrix: \n", id);
        print_matrix(array, 0, local_rows, local_cols);
        
        if (p > 1) {
            buffer = (int**)malloc(BLOCK_SIZE(p - 1, p, rows)) * sizeof(int*));
            buffStorage = (int*)malloc((BLOCK_SIZE(p - 1, p, rows)) * local_cols * sizeof(int));
            
            for (i = 0; i < (BLOCK_SIZE(p - 1, p, rows)); i++) {
                buffer[i] = &buffStorage[i * local_cols];
            }
            
            for(i = 1; i < p; i++) {
                MPI_Send(&prompt, 1, MPI_INT, i, PROMPT_MSG, comm);
                MPI_Recv(*buffer, local_rows * local_cols, MPI_INT, i, RESPONSE_MSG, comm, &status);
                
                print_matrix(buffer, 0, local_rows, local_cols);
            }
            
            printf("\n");
        }
    } else {
        MPI_Recv(&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm, &status);
        MPI_Send(*array, local_rows * local_cols, MPI_INT, 0, RESPONSE_MSG, comm);
    }
    
    return;
}


/*
Print Matrix
    parameters:
        array
        low_bound Customizes output for padded matrices.
        rows, cols; The boundaries for the loop.
*/
void print_matrix(int **array, int low_bound, int rows, int cols) {
    int i, j;
    
    for (i = low_bound; i < rows, i++) {
        for(j = low_bound; j < cols; j++) {
            printf("%d ", array[i][j]);
        }
        printf("\n");
    }
}