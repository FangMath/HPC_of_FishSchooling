/* This is the computational routine of mex file main routine FreeByFreeSmth.c
 * Copyright Fang Fang, 04/27/2015, Courant Inst.
 */
#include "mex.h"
#include "blas.h"
#include <omp.h>
#include <ctime>
#include <Eigen/Dense>
#include <iostream>
#define M_PI 3.14159265358979323846

using namespace Eigen;
using namespace std;
/* The computational routine */
void freebyfreesmth_single_routine(size_t nt, double *ztreal, double *ztimag, size_t ns, double *zsreal, double *zsimag, double *potreal, double *potimag, double *delta, double *distr)
{
    /* scalar values to use in dgemm */

    MatrixXd dx(nt,ns), dy(nt,ns), dz(nt,ns);
    VectorXd zztreal(nt), zztimag(nt);
    VectorXd ppotreal(nt), ppotimag(nt);
    VectorXd ddistr(ns);

    omp_set_num_threads(8);

    #pragma omp parallel 
    {
       int  pn = omp_get_num_threads();
    if (omp_get_thread_num() == 0) {
         printf("Number of threads = %d\n", omp_get_num_threads()); }
  printf("Thread %d is starting...\n",omp_get_thread_num());
    } 

#pragma omp parallel for
    for (int i = 0; i < nt; ++i){
    zztreal(i) = ztreal[i];
    zztimag(i) = ztimag[i];
    }

#pragma omp parallel for
    for (int i = 0; i < ns; ++i){
        ddistr(i) = distr[i];
    }

  /*
    #pragma omp parallel 
    {
       int  pn = omp_get_num_threads();
    if (omp_get_thread_num() == 0) {
         printf("Number of threads = %d\n", omp_get_num_threads()); }
  printf("Thread %d is starting...\n",omp_get_thread_num());
#pragma omp parallel for
    } 

  */

#pragma omp parallel for
    for (int j=0; j<nt; j++) {
#pragma omp parallel for
        for (int k=0; k<ns; k++) {
            dx(j,k)=(ztreal[j]-zsreal[k]);
            dy(j,k)=(ztimag[j]-zsimag[k]);
            dz(j,k)=dx(j,k)*dx(j,k)+dy(j,k)*dy(j,k)+delta[k]*delta[k];
            dx(j,k) = dx(j,k)/dz(j,k);
            dy(j,k) = dy(j,k)/dz(j,k);
            }
        }
    ppotreal = dx*ddistr;
    ppotimag = dy*ddistr;

#pragma omp parallel for
    for (int j=0; j<nt; j++) {
        potreal[j] = -ppotreal(j)/(2*ns*M_PI);
        potimag[j] = ppotimag(j)/(2*ns*M_PI);
    }
}
