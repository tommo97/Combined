#include "pgesv.hpp"
#include "array.hpp"


/**************************************************************/
void invrt(double* A, int N) {
    int *IPIV = new int[N + 1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N, &N, A, &N, IPIV, &INFO);
    dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

    delete [] IPIV;
    delete [] WORK;
}
/**************************************************************/
void inverse(Array <Array <double> > &Ain) {
    int sz = Ain.size();
    double *A = new double[sz * sz];
    int k = 0;
    for (int i = 0; i < sz; ++i)
        for (int j = 0; j < sz; ++j) {
            A[k] = Ain[i][j];
            k++;
        }

    //  Some memory problem occurring here...
    invrt(A, sz);

    k = 0;
    for (int i = 0; i < sz; ++i)
        for (int j = 0; j < sz; ++j) {
            Ain[i][j] = A[k];
            k++;
        }



    delete [] A;

}
/**************************************************************/
Array <double> pgesv(Array <Array <double> >&Ain, Array <double> &bin) {
    int sz = Ain.size();
    double **A = new double*[sz];

    double *rhs = new double[sz];
    double *gamma = new double[sz];

    for (int i = 0; i < sz; ++i) {
        A[i] = new double[sz];
        rhs[i] = bin[i];
        gamma[i] = 0.0;
    }

    for (int i = 0; i < sz; ++i)
        for (int j = 0; j < sz; ++j) {
            A[i][j] = Ain[i][j];
        }

    pgesv(A, rhs, sz, gamma);

    Array <double> gma(sz);

    for (int i = 0; i < sz; ++i)
        gma[i] = gamma[i];

    for (int i = 0; i < sz; ++i)
        delete [] A[i];

    delete [] A;
    delete [] gamma;
    delete [] rhs;

    return gma;


}
/**************************************************************/
#ifdef dgels
Array <double> dgels(Array < Array <double> > A, Array < double > X) {
    int M = A.size();
    if (M == 0)
        return 1;

    int N = A[M - 1].size();
    int nrhs = 1;
    int lda = M, ldb = M, info = 0, lwork = -1;
    double wkopt;
    double *work;

    double *a = arraytofortranarray(A);
    double *b = arraytofortranarray(X);

    char transp[] = "N"; // for "No transpose"

    dgels_(transp, &M, &N, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork,
            &info);

    lwork = (int) wkopt;
    work = new double[lwork];

    dgels_(transp, &M, &N, &nrhs, a, &lda, b, &ldb, work, &lwork,
            &info);

    if (info > 0) {
        printf("The diagonal element %i of the triangular factor ", info);
        printf("of A is zero, so that A does not have full rank;\n");
        printf("the least squares solution could not be computed.\n");
        exit(1);
    }
    delete [] work;
    delete [] a;
    delete [] b;
    return fortranarraytoarray(b, N);
}
#endif
/**************************************************************/
Array <double> fortranarraytoarray(double *in, int length){

    Array <double> out(length);
    for (int i = 0; i < length; i++)
        out[i] = in[i];
    return (out);
}



/**************************************************************/
double* arraytofortranarray(Array < Array <double> >  in) {
    int rows = in.size();
    int cols = in[0].size();
    double *out;
    int i, j;

    out = new double[rows * cols];
    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            out[i + j * cols] = in[i][j];

    return (out);
}


/**************************************************************/
double* arraytofortranarray(Array <double>  in) {
    double *out = new double[in.size()];
    for (int j = 0; j < in.size(); j++) out[j] = in[j];
    return (out);
}
/**************************************************************/
void pgesv(double **A, double *b, int n, double *x) {
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

    for (int i = 0; i < n; i++) x[i] = b[i];

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
/**************************************************************/
double* pgesv_ctof(double **in, int rows, int cols) {
    double *out;
    int i, j;

    out = new double[rows * cols];
    for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) out[i + j * cols] = in[i][j];
    return (out);
}
