/*
This file has my implementation of the LAPACK routine dgesv for C++.
This program solves the general set of double linear equations of the
form Ax = b, where A is a given n by n double matrix and b is a given
n element vector.  The n element vector x is returned.  The function
call is of the form

    void pgesv(double **A, double *b, int n, double *x)

    A: the left hand size n by n matrix
    b: the right hand side n element vector
    n: the dimension of A and b
    x: the n element vector for returned values

Scot Shaw
1 February 2000 */
#ifndef PGESV_HPP
#define PGESV_HPP
#include <math.h>
#include "array.hpp"
#define DOUBLE_PRECISION


void pgesv(double **A, double *b, int n, double *x);
double *pgesv_ctof(double **in, int rows, int cols);

#ifdef DOUBLE_PRECISION
extern "C" {
#ifdef dgels
    void dgels_( char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, double* work, int* lwork, int* info );
#endif
    //  Calc solution of general linear algebraic system
    void dgesv_(int *n, int *nrhs, double *a, int *lda,
            int *ipiv, double *b, int *ldb, int *info);

    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}
#else
extern "C" void sgesv_(int *n, int *nrhs, double *a, int *lda,
        int *ipiv, double *b, int *ldb, int *info);
#endif



void invrt(double* A, int N);

void inverse(Array <Array <double> > &Ain);

Array <double> pgesv(Array <Array <double> >&Ain, Array <double> &bin);


void dgels(double **A, double *b, int n, int m, double *x);

void pgesv(double **A, double *b, int n, double *x);

double* pgesv_ctof(double **in, int rows, int cols);

Array <double> fortranarraytoarray(double *in, int);

double* arraytofortranarray(Array < Array <double> >  in);

double* arraytofortranarray( Array <double>  in);
#ifdef dgels
Array <double> dgels(Array < Array <double> > A, Array < double > X);
#endif
#endif
