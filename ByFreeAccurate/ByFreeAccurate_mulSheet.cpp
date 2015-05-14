/*=========================================================
 * ByFreeAccurate_mulSheet.c - A mex routine that calculates the
 * complex velocity induced by free vortex sheet on itself,  
 *
 * Targets and sources are same points
 *
 * This is a MEX-file for MATLAB.
 * Copyright Fang Fang, 04/27/2015, Courant Inst.
 *=======================================================*/

#include "mex.h"
#include "blas.h"
#include <omp.h>


void ByFreeAccurate_mulSheet_routine(size_t nt, double *ztreal, double *ztimag, size_t ns, double *zetaf_r, double *zetaf_i, double *potreal, double *potimag, size_t Nw, double *L, double *gammaf);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ztreal, *ztimag, *zetaf_r, *zetaf_i, *potreal, *potimag, *L, *gammaf; /* pointers to input & output matrices*/
    size_t nt, ns, Nw;      /* vector dimensions */
    if(nrhs!=6) { mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Two inputs required.");
    }

    ztreal = mxGetPr(prhs[0]); /* target points real part */
    ztimag = mxGetPr(prhs[1]); /* target points imag part */
    zetaf_r = mxGetPr(prhs[2]); /* source points real part */
    zetaf_i = mxGetPr(prhs[3]); /* source points imag part */
    gammaf = mxGetPr(prhs[4]); /* circulation on free sheets */
    L = mxGetPr(prhs[5]); /* index of sheets */

    /* dimensions of input matrices */
    nt = mxGetM(prhs[0]);  
    ns = mxGetM(prhs[2]);  
    Nw = mxGetM(prhs[5]);  

    if (mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:first argument is not a vector!");
    }
    if (mxGetN(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:third argument is not a vector!");
    }
    if (mxGetN(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
 "MATLAB:sixth argument is not a vector!");
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
 "MATLAB:gammaf length does not match source length!");
    }

    /* create output matrix potreal and potimag */
    plhs[0] = mxCreateDoubleMatrix(nt, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nt, 1, mxREAL);
    potreal = mxGetPr(plhs[0]);
    potimag = mxGetPr(plhs[1]);

    /* Pass arguments to Fortran by reference */
    ByFreeAccurate_mulSheet_routine(nt, ztreal, ztimag, ns, zetaf_r, zetaf_i, potreal, potimag, Nw, L, gammaf);
}
