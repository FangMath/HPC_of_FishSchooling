/*=========================================================
 * FreeByFreeSmth.c - A mex routine that calculates the
 * complex velocity induced by free vortex sheet on itself,  
 * delta smoothing is applied.
 *
 * Targets and sources are same points
 *
 * This is a MEX-file for MATLAB.
 * Copyright Fang Fang, 04/27/2015, Courant Inst.
 *=======================================================*/

#include "mex.h"
#include "blas.h"
#include <omp.h>

void freebyfreesmth_single_routine(size_t nt, double *ztreal, double *ztimag, size_t ns, double *zsreal, double *zsimag, double *potreal, double *potimag, double *delta, double *distr);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ztreal, *ztimag, *zsreal, *zsimag, *potreal, *potimag, *delta, *distr; /* pointers to input & output matrices*/
    size_t nt, ns, n;      /* vector dimensions */
    if(nrhs!=6) { mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Two inputs required.");
    }

    ztreal = mxGetPr(prhs[0]); /* target points real part */
    ztimag = mxGetPr(prhs[1]); /* target points imag part */
    zsreal = mxGetPr(prhs[2]); /* source points real part */
    zsimag = mxGetPr(prhs[3]); /* source points imag part */
    delta = mxGetPr(prhs[4]); /* delta smoothing */
    distr = mxGetPr(prhs[5]); /* distribution on sources */

    /* dimensions of input matrices */
    nt = mxGetM(prhs[0]);  
    ns = mxGetM(prhs[2]);  

    if (mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:first argument is not a vector!");
    }
    if (mxGetN(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:third argument is not a vector!");
    }
    /* check input first dimensions */
    if ((mxGetM(prhs[1]) != nt) || mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:targets imag part length doesn't match real part!");
    }
    if ((mxGetM(prhs[3]) != ns) || mxGetN(prhs[3]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:sources imag part length doesn't match real part!");
    }
    if ((mxGetM(prhs[4]) != ns) || mxGetN(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:delta length does not match source length!");
    }
    if ((mxGetM(prhs[5]) != ns) || mxGetN(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:distr length does not match source length!");
    }

    /* create output matrix potreal and potimag */
    plhs[0] = mxCreateDoubleMatrix(nt, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nt, 1, mxREAL);
    potreal = mxGetPr(plhs[0]);
    potimag = mxGetPr(plhs[1]);

    /* Pass arguments to Fortran by reference */
    freebyfreesmth_single_routine(nt, ztreal, ztimag, ns, zsreal, zsimag, potreal, potimag, delta, distr);
}
