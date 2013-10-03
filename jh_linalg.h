#ifndef JH_LINALG
  #define JH_LINALG

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
//#include "JH_util.h"
#define Pi 3.1415926535897932385

/*needs to source libhoward and nrutil*/

void   AddAsgn_VV(double *, double *, int);
/*adds 2 arrays together
 * v1, returned with v1 + v2
 * v2, second vector
 * size of v1 and v2
 */
void   Asgn_VV(double *, double *, int);
/* copies an array to another array
 * v1 is null, returned with v2
 * v2, vector to be copied to v1
 * size of v1 and v2
 */
void   gaussj(double **, int, double **, int);
/*calculates the inverse of a matrix using gaussian elimination
 * m1 input, output Inverse(m1)
 * rank of m1
 * v1,normalization multiplies for each row
 * size of v1, usely 1
 */
void   JMat_Inv_2x2(double **, double **);
/*calculates the inverse of a matrix of rank 2
 * m1^-1 to be returned
 * m2, input
 */
double   l2Norm_V(double *, int);
/*returns the sqrt of (v1.v1),
 * v1, array
 * n, size of v1
 */
void   lubksb(double **, int, int *, double *);
/* calculates solution to x to Ax=b
 * A, matrix
 * n, rank and sizeo of A and b
 * v1, contains record of row swapping
 * b input, output x
 */
void   ludcmp(double **, int, int *, double *);
/*decomposes a matrix into LU part
 * m1 on in, returns LU out
 * n, rank of m1
 * v1, contains vector to swap rows in m1 to place max value in diagonol
 * p, scalar 1 or -1 containing the permutations for matrix
 */
void   Matrix_Inverse(double **,  double **, int);
/* calculates the inverse of a matrix using LU decomposition
 * m1, returned with Inverse of m2
 * m2, input matrix
 * n, rank of m1 and m2
 */
void   Matrix_Inverse2(double **,  double **, int);
/* calculates the inverse of a matrix using LU decomposition
 * m1, returned with Inverse of m2
 * m2, input matrix
 * n, rank of m1 and m2
 */
void   Matrix_Matrix(double **,  double **,  double **, int);
/*multiplies 2 matrices together and returns in a 3rd, m3 = m1.m2
 * m1 is null, returned with m2.m3
 * m2, first matrix
 * m3, 2nd matrix in
 * n, rank of m1, m2, m3
 */
void   JH_inp1_mat_mat(double **,  double **, int);
/*multiplies 2 matrices together and returns 1st mat
 * m1, first matrix, also returned matrix
 * m2, 2nd matrix in
 * n, rank of m1, m2, m3
 */
void   JH_inp2_mat_mat(double **,  double **, int);
/*multiplies 2 matrices together and returns 1st mat
 * m1, first matrix
 * m2, 2nd matrix, also returned matrix
 * n, rank of m1, m2, m3
 */
void   Matrix_Vec(double *,  double **,  double *, int, int);
/* multiplies a matrix times a vector, v1 = m1.v2
 * v1
 * m1, matrix
 * v2, input vector
 * nr, number of rows of m1
 * nc, number of columns of m1 ans size of v1 and v2
 */
double * Mul_MV(double **, double *, int, int);
/* multiplies a matrix times a vector, return = m1.v1
 * m1, matrix
 * v1, input vector
 * nr, number of rows of m1
 * nc, number of columns of m1 ans size of v1
 */
double * Mul_SV(double, double *, int);
/*returns the an array multiplied by an scalar, return = s*v1
 * s, scalar
 * v1, array
 * n,size of v1
 */
double   Mul_VV(double *, double *, int);
/*return the value of return = v1.v2
 * v1, array
 * v2, array
 * n, size of v1 and v2
 */
double * MulAsgn_VS(double *, double, int);
/*returns and reasigns v1 the value of (return,v1) = s*v1
 * v1, array, returned with s*v1
 * s, scalar value
 * n, size of v1
 */
void   Newton_Solver(double *,  double **,  double *, int);
/*calculates the solution x to Ax=b using LU decomposition
 * x, solution
 * A, input matrix
 * b, nonhomogenious part of Ax=b
 * rank and size of A, x, b
 */
void   Set_VV(double *, double *, int);
/*copies v2 to v1
 * v1
 * v2
 * size of v1, v2
 */
void   Sub_VV(double *, double *, double *, int);
/* subtracts 2 arrays, v1 = v2 - v3
 * v1
 * v2
 * v3
 * size of v1,v2,v3
 */
double * SubAsgn_VV(double *, double *, int);
/* returns and reasigns v1 = v1-v2
 * v1, array, returned with v1 - v2
 * v2
 * size of v1 and v2
 */
double   Vec_Vec(double *,  double *, int);
/*returns the values of v1.v2
 * v1
 * v2
 * size of v1 and v2
 */
#endif

