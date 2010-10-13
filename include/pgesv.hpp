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
#define DOUBLE_PRECISION


void pgesv(double **A, double *b, int n, double *x);
double *pgesv_ctof(double **in, int rows, int cols);

#ifdef DOUBLE_PRECISION
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda,
                    int *ipiv, double *b, int *ldb, int *info);
#else
extern "C" void sgesv_(int *n, int *nrhs, double *a, int *lda,
                    int *ipiv, double *b, int *ldb, int *info);
#endif


Array <double> pgesv(Array <Array <double> >&Ain, Array <double> &bin)
{
    int sz = Ain.size();
    double **A = new double*[sz];

    double *rhs = new double[sz];
    double *gamma = new double[sz];

    for (int i = 0; i < sz; ++i)
    {
        A[i] = new double[sz];
        rhs[i] = bin[i];
        gamma[i] = 0.0;
    }

    for (int i = 0; i < sz; ++i)
        for (int j = 0; j < sz; ++j){
            A[i][j] = Ain[i][j];
        }

    pgesv(A, rhs, sz, gamma);

    Array <double> gma(sz);

    for (int i = 0; i < sz; ++i)
        gma[i] = gamma[i];

        for (int i = 0; i < sz; ++i)
                delete [] A[i];

        delete [] A;
        delete gamma;
        delete rhs;

    return gma;


}

void pgesv(double **A, double *b, int n, double *x)
{
int nrhs, lda, ldb, info;
double *a;
int *ipiv;

nrhs = 1; /* The Fortran routine that I am calling will solve a
            system of equations with multiple right hand sides.
            nrhs gives the number of right hand sides, and then b
            must be a matrix of dimension n by nrhs. */

lda = n; // The leading dimension of A
a = pgesv_ctof(A, lda, n); // Convert A to a Fortran style matrix

ipiv = new int[n]; // Allocate memory to store pivot information

ldb = n;

/* The Fortran routine replaces the input information in b with the
    results.  Since we are interested in the results, we put the
    initial conditions in the output array and use that in the
    call. */

for (int i=0; i<n; i++) x[i] = b[i];

// Now the function call
#ifdef DOUBLE_PRECISION
dgesv_(&n, &nrhs, a, &lda, ipiv, x, &ldb, &info);
#else
sgesv_(&n, &nrhs, a, &lda, ipiv, x, &ldb, &info);
#endif
// Clean up the memory before returning to the calling program

delete a;
delete ipiv;
}


double* pgesv_ctof(double **in, int rows, int cols)
{
double *out;
int i, j;

out = new double[rows*cols];
for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*cols] = in[i][j];
return(out);
}
#endif