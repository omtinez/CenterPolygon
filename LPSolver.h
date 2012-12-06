/*
 * Author: Oscar Martinez
 * Copyright (c) 2012 All rights reserved.
 */

#ifndef LPSOLVER_H_
#define LPSOLVER_H_

typedef struct Vector Vector;
struct Vector {
    unsigned int  n;
    float*			values;
};

typedef struct Matrix_2D Matrix_2D;
struct Matrix_2D {
    unsigned int	rows;
    unsigned int	cols;
    float*			values;
};


Vector* Vector_new(unsigned int);
void Vector_destroy(Vector*);

Matrix_2D* Matrix_2D_new(int, int);
void Matrix_2D_destroy(Matrix_2D*);

float Vector_get(Vector*, int);
void Vector_set_val(Vector*, int, float);
void Vector_set_all(Vector*, double, ... );
void Vector_print(Vector* vector);

float Matrix_2D_get_cell(Matrix_2D*, int, int);
void Matrix_2D_set_cell(Matrix_2D*, int, int, float);
void Matrix_2D_set_row(Matrix_2D*, int, double, ... );
void Matrix_2D_set_all(Matrix_2D*, double, ... );
void Matrix_2D_scale_row(Matrix_2D*, unsigned int, float);
void Matrix_2D_print(Matrix_2D*);

Vector* LPSolver_simplex_solve(Matrix_2D*);

void LPSolver_test_vector(void);
void LPSolver_test_simplex_1(void);
void LPSolver_test_simplex_2(void);
void LPSolver_test_simplex_3(void);



#endif /* LPSOLVER_H_ */