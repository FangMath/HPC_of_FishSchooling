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

void freebyfreesmth_mul_routine(size_t nt, double *ztreal, double *ztimag, double *Lt,  size_t ns, double *zsreal, double *zsimag, double *potreal, double *potimag, double *delta, double *distr, double *Ls, size_t Nw);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ztreal, *ztimag, *Lt, *zsreal, *zsimag, *delta, *distr, *Ls, *potreal, *potimag; /* pointers to input & output matrices*/
    size_t nt, ns, Nw;      /* vector dimensions */
    if(nrhs!=8) { mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Eight inputs required.");
    }

    ztreal = mxGetPr(prhs[0]); /* target points real part */
    ztimag = mxGetPr(prhs[1]); /* target points imag part */
    Lt = mxGetPr(prhs[2]); /* index of sheets */
    zsreal = mxGetPr(prhs[3]); /* source points real part */
    zsimag = mxGetPr(prhs[4]); /* source points imag part */
    delta = mxGetPr(prhs[5]); /* delta smoothing */
    distr = mxGetPr(prhs[6]); /* distribution on sources */
    Ls = mxGetPr(prhs[7]); /* index of sources */

    /* dimensions of input matrices */
    nt = mxGetM(prhs[0]);  
    Nw = mxGetM(prhs[2]);  
    ns = mxGetM(prhs[3]);  

    if (mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:first argument is not a vector!");
    }
    if (mxGetN(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:third argument is not a vector!");
    }
    if (mxGetN(prhs[3]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:fourth argument is not a vector!");
    }

    /* check input first dimensions */
    if ((mxGetM(prhs[1]) != nt) || mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:targets imag part length doesn't match real part!");
    }
    if ((mxGetM(prhs[4]) != ns) || mxGetN(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:sources imag part length doesn't match real part!");
    }
    if ((mxGetM(prhs[5]) != ns) || mxGetN(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:delta length does not match source length!");
    }
    if ((mxGetM(prhs[6]) != ns) || mxGetN(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:distr length does not match source length!");
    }
    if ((mxGetM(prhs[7]) != Nw) || mxGetN(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:Ls length does not match Lt length!");
    }

    /* create output matrix potreal and potimag */
    plhs[0] = mxCreateDoubleMatrix(nt, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nt, 1, mxREAL);
    potreal = mxGetPr(plhs[0]);
    potimag = mxGetPr(plhs[1]);

    /* Pass arguments to Fortran by reference */
    freebyfreesmth_mul_routine(nt, ztreal, ztimag, Lt, ns, zsreal, zsimag, potreal, potimag, delta, distr, Ls, Nw);
}
