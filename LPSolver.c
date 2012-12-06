/*
 * Author: Oscar Martinez
 * Copyright (c) 2012 All rights reserved.
 */

// library includes
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <float.h>

// module includes
#include "LPSolver.h"


// vector constructor
Vector* Vector_new(unsigned int n) {
    // memory allocation
    Vector* this = malloc(sizeof(Vector));
    this->values = malloc(sizeof(float)*(size_t)n);
    // populate object
    this->n = n;
    // initialize vector
    int i;
    for(i=0;i<n;i++) this->values[i]=0;

    return this;
}

// vector destructor
void Vector_destroy(Vector* vector) {
    // memory de-allocation
  if (vector == NULL) return;
	free (vector->values);
	free (vector);
}

// function used to retrieve a value (n) of a vector
float Vector_get(Vector* vector, int n) {
    return vector->values[n];
}

// function used to set a value (n) of a vector
void Vector_set_val(Vector* vector, int n, float val) {
	vector->values[n] = val;
}

// convenience function that sets an entire vector in one sweep
void Vector_set_all(Vector* vector, double val, ... ) {
    va_list args;
    va_start ( args, val );

    int i;
//    float tmp = 0;
    for (i=0; i<vector->n; i++) {
    	Vector_set_val(vector, i, (float) val);
        val = va_arg(args, double);
    }
//    printf("total: %f\n",tmp);
   	va_end( args );
}

// multiply the whole vector by a scalar
void Vector_scale(Vector *vector, float k) {
    int i;
    for (i=0; i<vector->n; i++)
    	vector->values[i] *= k;
//        (vector->values)[sizeof(float)*i] *= k;
}


// vector toString()
void Vector_print(Vector* vector) {
	int i;
	for (i=0; i<vector->n; i++) {
		if (i==0) printf("[%f",Vector_get(vector,i));
		else printf(",%f",Vector_get(vector,i));
	}	printf("]\n");
}

// matrix constructor
Matrix_2D* Matrix_2D_new(int rows, int cols) {
    // memory allocation
    Matrix_2D* this = malloc(sizeof(Matrix_2D));
    this->values = malloc(sizeof(float) * (size_t) rows * (size_t) cols);
    // populate object
    this->rows = rows;
    this->cols = cols;
    // initialize matrix
    int i,j;
    for(i = 0; i < rows; i++) for(j = 0; j < cols; j++) (this->values)[rows * j + i] = 0;

    return this;
}

// matrix destructor
void Matrix_2D_destroy(Matrix_2D *matrix) {
    // memory de-allocation
	free (matrix->values);
	free (matrix);
}

// function used to retrieve a value (m,n) of a 2d matrix
float Matrix_2D_get_cell(Matrix_2D *matrix, int row, int col) {
    return matrix->values[row * matrix->cols + col];
}

// function used to retreive as a vector a whole row of a 2d matrix
Vector* Matrix_2D_get_row(Matrix_2D *matrix, int row) {
    Vector* vector = Vector_new(matrix->cols);

    // get the values of the row from the matrix and put them in the vector
    int i;
    for(i = 0; i < matrix->cols; i++)
    	Vector_set_val(vector, i, Matrix_2D_get_cell(matrix, row, i));

    return vector;
}

// function used to retreive as a vector a whole column of a 2d matrix
Vector* Matrix_2D_get_col(Matrix_2D *matrix, int col) {
    Vector* vector = Vector_new(matrix->rows);

    // get the values of the row from the matrix and put them in the vector
    int i;
    for(i = 0 ; i < matrix->rows; i++) Vector_set_val(vector, i, Matrix_2D_get_cell(matrix, i, col));

    return vector;
}

// function used to set a value (m,n) of a 2d matrix
void Matrix_2D_set_cell(Matrix_2D *matrix, int row, int col, float val) {
	matrix->values[row * matrix->cols + col] = val;
}

// convenience function that sets entire rows of a matrix in one sweep
void Matrix_2D_set_row(Matrix_2D *matrix, int row, double val, ... ) {
    va_list args;
    va_start ( args, val );

    int i;
    for (i = 0; i < matrix->cols; i++) {
        Matrix_2D_set_cell(matrix, row, i, (float) val);
        val = va_arg(args, double);
    }

    va_end(args);
}

// convenience function that sets entire rows of a matrix in one sweep
void Matrix_2D_set_all(Matrix_2D *matrix, double val, ... ) {
    va_list args;
    va_start ( args, val );

    int i,j;
    for(i=0;i<matrix->rows;i++) for(j=0;j<matrix->cols;j++) {
    	Matrix_2D_set_cell(matrix, i, j, (float) val);
    	val = va_arg (args, double);
    }

    va_end(args);
}

// convenience function that sets entire rows of a matrix in one sweep
//void Matrix_2D_set_vector(Matrix_2D *matrix, unsigned int m, Vector* vector ) {
//    int i;
//    for (i=0; i<vector->n; i++)
//        Matrix_2D_set_cell(matrix,m,i,Vector_get(vector,i));
//}

// multiply an entire row of the matrix by a scalar
void Matrix_2D_scale_row(Matrix_2D *matrix, unsigned int row, float k) {
    int i;
    for (i=0; i<matrix->cols; i++)
        (matrix->values)[row * matrix->cols + i] *= k;
}

// sum row m2 * scalar to row m1
void Matrix_2D_sum_rows(Matrix_2D *matrix, unsigned int row1, unsigned int row2, float k) {
    int i;
    for (i=0; i<matrix->cols; i++)
        matrix->values[row1 * matrix->cols + i] +=
            k * Matrix_2D_get_cell(matrix, row2, i);
}

// matrix toString()
void Matrix_2D_print(Matrix_2D* matrix) {
	int i,j;
	printf("[");
	for (i=0; i<matrix->rows; i++) {
		for (j=0; j<matrix->cols; j++) {
			if (j==0) printf("[%f",Matrix_2D_get_cell(matrix,i,j));
			else printf(",%f",Matrix_2D_get_cell(matrix,i,j));
		}	if (i!=matrix->rows-1) printf("]\n");
	}	printf("]]\n");
}

Vector* LPSolver_simplex_solve(Matrix_2D* problem_data) {
	// declare all the working variables for this function
	int i, j, pivot_row = 0,pivot_col = 0, unbounded_flag;
	float tmp, min_value;

	// create a matrix based on the problem_data input matrix but include slack variables
	Matrix_2D *tableau = Matrix_2D_new(problem_data->rows, problem_data->cols + problem_data->rows - 1);
	for(i = 0; i < tableau->rows; i++) {
		for(j = 0; j < tableau->cols; j++) {

			// when the row is the cost function
			if (i == tableau->rows - 1) {
				Matrix_2D_set_cell(tableau, i, j,
					Matrix_2D_get_cell(problem_data, i, j)*-1);

			// when the column is a normal variable
			} else if(j < problem_data->cols - 1) {
				Matrix_2D_set_cell(tableau, i, j,
					Matrix_2D_get_cell(problem_data, i, j));

			// when the column is a slack variable
			} else if(j < tableau->cols - 1) {
				if (j == i + problem_data->cols - 1)
					Matrix_2D_set_cell(tableau, i, j, 1.0);
				else
					Matrix_2D_set_cell(tableau, i, j, 0.0);

			// when the column is the independent term
			} else if(j == tableau->cols - 1) {
				Matrix_2D_set_cell(tableau, i, tableau->cols - 1,
					Matrix_2D_get_cell(problem_data,i,problem_data->cols - 1));
			}
		}
	}
    
//    Matrix_2D_print(tableau);

    // label for the beginning of the simplex method
//	int count=0;
    simplex_method:
//    printf("\nTableau, iteration %d:\n", ++count);
//    Matrix_2D_print(tableau);
    
    // search for unbounded solution case
    unbounded_flag = 0;
    for (i = 0; i < tableau->cols; i++)
    	if (Matrix_2D_get_cell(tableau, tableau->rows - 1, i) > 0)
    		for (j = 0; j < tableau->rows; j++)
    			if (Matrix_2D_get_cell(tableau, j, i) >= 0)
    				unbounded_flag++;
    if (unbounded_flag == tableau->rows - 1) return NULL;

    // find smallest coefficient in cost function to set pivot column
    min_value = DBL_MAX;
    for(i = 0; i < tableau->cols; i++) {
        tmp = Matrix_2D_get_cell(tableau, tableau->rows - 1, i);
        if (tmp < min_value) {
            min_value = tmp;
            pivot_col = i;
        }
    }

    // if all coefficients are positive then this is the final tableau
    if((float) min_value >= 0) {
//		Matrix_2D_print(tableau);

    	// if one of the artificial or slack variables are still positive then it's infeasible solution
    	for (i = tableau->rows - 1; i < tableau->cols - 1; i++)
    		for (j = 0; j < tableau->rows; j++)
    			if (Matrix_2D_get_cell(tableau, i, j) < 0.0f) return NULL;

    	//

    	// find the solution for the basic variables from the tableau
    	Vector* solution = Vector_new(tableau->cols - 1);
    	for (i = 0; i < tableau->rows - 1; i++)
    		for (j = 0; j < tableau->cols - 2; j++)
    			if (Matrix_2D_get_cell(tableau, i, j) == 1.0f)
    				Vector_set_val(solution, j, Matrix_2D_get_cell(tableau, i, tableau->cols - 1));

    	// cleanup and return
    	Matrix_2D_destroy(tableau);
    	return solution;
    }
    
    // otherwise find minimum positive ratio between the variable from the pivot column
    // and the independent term to determine the pivot row
    min_value = DBL_MAX;
    for(i = 0; i < tableau->rows - 1; i++) {
		tmp = Matrix_2D_get_cell(tableau, i, pivot_col);
		if ((float)tmp >= 0) {
			tmp = Matrix_2D_get_cell(tableau, i, tableau->cols - 1) / tmp;
			if (tmp < min_value) {
				min_value = tmp;
				pivot_row = i;
			}
        }
    }
//    printf("pivot col: %d, pivot row: %d (ratio: %f)\n",pivot_col,pivot_row,min_value);
    

    // make the coefficient for the pivot element 1
    tmp = 1 / Matrix_2D_get_cell(tableau, pivot_row, pivot_col);
    Matrix_2D_scale_row(tableau, pivot_row, tmp);
    
    // perform the simplex pivot operation by making all other elements in pivot column equal to zero
    for(i = 0; i < tableau->rows; i++) {
        if (i != pivot_row) {
        	Matrix_2D_sum_rows(tableau, i, pivot_row,
                Matrix_2D_get_cell(tableau, i, pivot_col) * -1);
        }
    }
    
    // continue performing the simplex method until it returns
    goto simplex_method;
    return 0; // never reaches this, but compiler cries otherwise
}



// UNIT TESTING
void LPSolver_test_vector(void) {
	Vector* v = Vector_new(10);
	Vector_set_all(v, 1.1, 2.2, 3.3f, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9f, 10);
	Vector_set_val(v, 9, 0.1);

//	printf("%s[TEST] ", timestamp());
	Vector_print(v); // [1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,0.1]
	Vector_destroy(v);
}

void LPSolver_test_simplex_1(void) {
	Matrix_2D* m = Matrix_2D_new(3,4);
	Matrix_2D_set_all(m,	0.5, 2.0, 1.0, 24.0,
							1.0, 2.0, 4.0, 60.0,
							6.0, 14.0, 13.0, 0.0 );

	Vector* v = LPSolver_simplex_solve(m);
//	printf("%s[TEST] Results: ", timestamp());
	Vector_print(v);

	Vector_destroy(v);
	Matrix_2D_destroy(m);
}

void LPSolver_test_simplex_2(void) {
	Matrix_2D* m = Matrix_2D_new(5,3);
	Matrix_2D_set_all(m,	8.0, 4.0, 3440.0,
							6.0, 8.0, 2880.0,
							4.0, 6.0, 2760.0,
							1.0, 0.0, 420.0,
							14.0, 16.0, 0.0 );

	Vector* v = LPSolver_simplex_solve(m);
//	printf("%s[TEST] Results: ", timestamp());
	Vector_print(v);

	Vector_destroy(v);
	Matrix_2D_destroy(m);
}

void LPSolver_test_simplex_3(void) {
	Matrix_2D* m = Matrix_2D_new(4,4);
	Matrix_2D_set_all(m,	1.0, 2.0, 2.0, 20.0,
							2.0, 1.0, 2.0, 20.0,
							2.0, 2.0, 1.0, 20.0,
							10.0, 12.0, 12.0, 0.0);

	Vector* v = LPSolver_simplex_solve(m);
//	printf("%s[TEST] Results: ", timestamp());
	Vector_print(v);

	Vector_destroy(v);
	Matrix_2D_destroy(m);
}