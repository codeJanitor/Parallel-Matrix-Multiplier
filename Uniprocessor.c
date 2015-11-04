/**
 * @author Caleb Mays and Nathan Brown
 * @date 16 Nov 2010
 * Description: Parallel - project 4, part 1. Matrix Multiplication.
 *              Provides a uniprocessor implementation of matrix multiplication.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <mpi.h>

#define BUFFER_SIZE 10

int main (int argc, char *argv[]) {
	// columns of A must equal rows of B
	// C will have rowsA x colsB
	int *rowsA, *colsA, *rowsB, *colsB, *rowsC, *colsC;
	int **a;
	int *aStorage;
	int **b;
	int	*bStorage;
	int **c;
	int *cStorage;
	int status;
	int i, j, k;
    
	if (argc != 2){
             printf("not correct set of arguments");
             return 0;
    }
    
	int id;
	int p;
	double elapsed_time;
	
	MPI_Init(&argc, &argv);	
	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
    char numBuffer[BUFFER_SIZE];
    
	// get the memory I need.
	rowsA = (int*)malloc(sizeof(int));
	colsA = (int*)malloc(sizeof(int));
	rowsB = (int*)malloc(sizeof(int));
	colsB = (int*)malloc(sizeof(int));	
	rowsC = (int*)malloc(sizeof(int));
	colsC = (int*)malloc(sizeof(int));
	
	if (id == (p-1)){
		// this is going to be all square. also, a power of two.
		*colsB = *colsC= *colsA = *rowsB = *rowsA = *rowsC = atoi(argv[1]);
	}


	a = (int**)malloc((*rowsA) * sizeof(int*));
	b = (int**)malloc((*rowsB) * sizeof(int*));
	c = (int**)malloc((*rowsC) * sizeof(int*));
	
	aStorage = (int*)malloc((*rowsA) * (*colsA) * sizeof(int));
	bStorage = (int*)malloc((*rowsB) * (*colsB) * sizeof(int));
	cStorage = (int*)malloc((*rowsC) * (*colsC) * sizeof(int));
	
	
	for (i = 0; i < (*rowsA); i++) {
		a[i] = &aStorage[i * (*rowsA)];
	}
	for (i = 0; i < (*rowsB); i++) {
		b[i] = &bStorage[i * (*rowsB)];
	}
	for (i = 0; i < (*rowsC); i++) {
		c[i] = &cStorage[i * (*rowsC)];
	}	
	
	for (i = 0; i < (*rowsA); i++) {
		for(j = 0; j < (*colsA); j++){
			a[i][j] = rand() % 11;
			b[i][j] = a[i][j];
			c[i][j] = 0;
		}
	}
	
	printf("\n");
	printf("\n");
	elapsed_time = -MPI_Wtime();
	for(i = 0; i < *rowsA; i++){
		for(j = 0; j < *colsB; j++){
			c[i][j] = 0;
			for (k = 0; k < *colsA; k++){
				c[i][j] += a[i][k] * b[k][j];				 
			}
		}
	}
    
	elapsed_time += MPI_Wtime();
    
    // Processor 0 writes results to the file.
	if (id == 0) {
		FILE *file;
		file = fopen("myAnswers.csv", "a+");
		fprintf(file, "%d, %d, %10.6f\n", p, *rowsA, elapsed_time);
		fclose(file);
	}
    
    // Free memory and exit.
	free(a);
	free(aStorage);
	free(b);
	free(bStorage);
	free(c);
	free(cStorage);
	free(rowsA);
	free(rowsB);
	free(rowsC);
	free(colsA);
	free(colsB);
	free(colsC);
	MPI_Finalize();
    
	return 0;
}
