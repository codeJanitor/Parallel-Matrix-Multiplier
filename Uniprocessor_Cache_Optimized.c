/**
 * @author Caleb Mays and Nathan Brown
 * @date 16 Nov 2010
 * Description: Parallel - project 4, part 2. Matrix Multiplication.
 *              Provides a uniprocessor implementation of matrix multiplication,
 *              but recursively splits the matrices to optimize caching performance.
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <sys/stat.h>
#include <mpi.h>


#define BUFFER_SIZE 10


void MatrixMultiply(int rowA, int **cMatrix, int **aMatrix, int **bMatrix, int crow, int ccol, int arow, int acol, int brow, int bcol);

int main (int argc, char *argv[]) {
    
	// columns of A must equal rows of B
	// C will have rowsA x colsB
	int *rowsA, *colsA, *rowsB, *colsB, *rowsC, *colsC;
	int     **a;     
	int     *aStorage;
	int 	**b;
	int		*bStorage;
	int 	**c;
	int 	*cStorage;
	int      status;
	int i, j, k;
    
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
    
	// get the memory I need.
	rowsA = (int*)malloc(sizeof(int));
	colsA = (int*)malloc(sizeof(int));
	rowsB = (int*)malloc(sizeof(int));
	colsB = (int*)malloc(sizeof(int));	
	rowsC = (int*)malloc(sizeof(int));
	colsC = (int*)malloc(sizeof(int));
	
	if (id == (p-1)){
		inputFile = fopen (inputFileName, "r");
        
		if (inputFile == NULL){
				   printf("No file pointer!\n");
				   return 0;
		}
		else{
			// this is going to be all square. also, a power of two.
			*colsB = *colsC = *colsA = *rowsB = *rowsA = *rowsC = atoi(fgets(numBuffer, BUFFER_SIZE, inputFile));
		}
	}

	// all of the memory allocation business and addressing business.
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
	

	// reading in and assigning the array.
	// I've simply decided that we're going to 
	// multiply a matrix against itself, so that's why
	// b will equal a.
	// also, i have to 0 out c.
	for (i = 0; i < (*rowsA); i++) {
		for(j = 0; j < (*colsA); j++){
              // all i really want to do is read in from a file,
              // and i don't care if i douplicate the other array
              // because i'm just simply running tests on sizes,
              // not for different numbers like with sorting or something.
			a[i][j] = fgetc(inputFile) - 48;
			b[i][j] = a[i][j];
			c[i][j] = 0;
		}
		fgetc(inputFile);
	}

	// time it, do it, stop timing.
	elapsed_time = -MPI_Wtime();
	MatrixMultiply(*rowsA, c, a, b, 0, 0, 0, 0, 0, 0);
	elapsed_time += MPI_Wtime();
    
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

/*
void MatrixMultiply(int l, int m, int n, int **cMatrix, int **aMatrix, int **bMatrix, int Rowc,
     int Colc, int Rowa, int Cola, int Rowb, int Colb){
    */
void MatrixMultiply(int rowA, int **cMatrix, int **aMatrix, int **bMatrix,
                    int crow, int ccol, int arow, int acol, int brow, int bcol) {
    
	int i, j, k;
	int colA, rowB, colB, rowC, colC;
	colA = rowB = rowB = colB = rowC = colC = rowA;
	int half[2];
	
	if (brow * bcol > 256) {
		// split this thing in half.
		half[0] = 0;
		half[1] = rowA/2;
        
		for(i = 0; i < 2; i++){
            for (j = 0; j < 2; j++){
                for (k = 0; k < 2; k++){
                    MatrixMultiply((rowA-rowA/2), cMatrix, aMatrix, bMatrix, crow + half[i], rowA + half[i], acol + half[k], brow + half[k], bcol + half[j]);
                }
			}
        }
	}
	
	
	else {
    	for(i = 0; i < rowA; i++) {
    		for(j = 0; j < colB; j++) {
    			cMatrix[i][j] = 0;
                
    			for (k = 0; k < colA; k++) {
                    cMatrix[i][j] += aMatrix[i][k] * bMatrix[k][j];
    			}
    		}
    	}
	}
    
	return;
}