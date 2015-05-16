/* This is the computational routine of mex file main routine FreeByFreeSmth.c
 * Copyright Fang Fang, 04/27/2015, Courant Inst.
 */
#include "mex.h"
#include "blas.h"
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <Eigen/Dense>
#include <Fish>
#include <iostream>
#include <vector>
#define M_PI 3.14159265358979323846

using namespace Eigen;
using namespace Fish;
using namespace std;

void freeSmth(VectorXd zetaf_r, VectorXd zetaf_i, VectorXd source_r, VectorXd source_i, VectorXd distr, VectorXd delta, VectorXd& pot_r, VectorXd& pot_i);

void freebyfreesmth_mul_routine(size_t nt, double *ztreal, double *ztimag, double *Lt, size_t ns, double *zsreal, double *zsimag, double *potreal, double *potimag, double *delta, double *distr, double *Ls, size_t Nw)
{
    /* scalar values to use in dgemm */

    VectorXd zzetaf_r(nt), zzetaf_i(nt);
    VectorXd ssource_r(ns), ssource_i(ns);
    VectorXd ddistr(ns), ddelta(ns);
    
/*
#pragma omp parallel 
{
int  pn = omp_get_num_threads();
if (omp_get_thread_num() == 0) {
printf("Number of threads = %d\n", omp_get_num_threads()); }
printf("Thread %d is starting...\n",omp_get_thread_num());
} 
*/

    /* read data from array to eigen vectors */
    for (int i = 0; i < nt; ++i){
        zzetaf_r(i) = ztreal[i];
        zzetaf_i(i) = ztimag[i];
    }

for (int i = 0; i < ns; ++i){
    ssource_r(i) = zsreal[i];
    ssource_i(i) = zsimag[i];
    ddistr(i) = distr[i];
    ddelta(i) = delta[i];
}

vector<Wing> W(Nw);
VectorXd pot_r = VectorXd::Zero(nt,1), pot_i = VectorXd::Zero(nt,1);

/* Computing */
cout << "max threads=" << omp_get_max_threads() << endl;
int Np = Nw;
omp_set_num_threads(Np);
#pragma omp parallel for 
for (int ib = 0; ib < Nw; ++ib){
    /*
printf("Number of threads = %d\n", omp_get_num_threads()); 
      printf("Thread %d is starting...\n",omp_get_thread_num());
      */

    /* read free sheets from array to c++ class */
    if (ib == 0){
        W[0].zetaf_r = zzetaf_r.head(Lt[0]);
        W[0].zetaf_i = zzetaf_i.head(Lt[0]);

        W[0].source_r = ssource_r.head(Ls[0]);
        W[0].source_i = ssource_i.head(Ls[0]);
        W[0].distr = ddistr.head(Ls[0]);
        W[0].delta = ddelta.head(Ls[0]);
        W[0].pot_r = VectorXd::Zero(Ls[0],1); 
        W[0].pot_i = VectorXd::Zero(Ls[0],1);
    }
    else
    {
        W[ib].zetaf_r = zzetaf_r.segment(Lt[ib-1], Lt[ib]-Lt[ib-1]);
        W[ib].zetaf_i = zzetaf_i.segment(Lt[ib-1], Lt[ib]-Lt[ib-1]);
        W[ib].source_r = ssource_r.segment(Ls[ib-1], Ls[ib]-Ls[ib-1]);
        W[ib].source_i = ssource_i.segment(Ls[ib-1], Ls[ib]-Ls[ib-1]);
        W[ib].distr = ddistr.segment(Ls[ib-1], Ls[ib]-Ls[ib-1]);
        W[ib].delta = ddelta.segment(Ls[ib-1], Ls[ib]-Ls[ib-1]);
        W[ib].pot_r = VectorXd::Zero(Ls[ib]-Ls[ib-1],1); 
        W[ib].pot_i = VectorXd::Zero(Ls[ib]-Ls[ib-1],1);
    }

    freeSmth(W[ib].zetaf_r, W[ib].zetaf_i, W[ib].source_r, W[ib].source_i, W[ib].distr, W[ib].delta, W[ib].pot_r, W[ib].pot_i);

    /* Convert output from eigen vector to array */
    if (ib == 0){
        for (int j=0; j<Lt[0]; j++) {
            potreal[j] = W[0].pot_r(j);
            potimag[j] = W[0].pot_i(j);
        }
    }
    else{
        for (int j=0; j<Lt[ib]-Lt[ib-1]; j++) {
            int idx = Lt[ib-1]+j;
            potreal[idx] = W[ib].pot_r(j);
            potimag[idx] = W[ib].pot_i(j);
        }
    }
}

}



/******************************* 
 * Implementation of freeSmth  *
 *******************************/
void freeSmth(VectorXd zetaf_r, VectorXd zetaf_i, VectorXd source_r, VectorXd source_i, VectorXd distr, VectorXd delta, VectorXd& pot_r, VectorXd& pot_i){

    /*
       timeval time1, time2;
       double ttl1, ttl2, diffl12;
       gettimeofday(&time1, NULL); ttl1 = (time1.tv_sec * 1000000 + time1.tv_usec);
       */
    int nt = zetaf_r.size(), ns = source_r.size();

    MatrixXd dx(nt,ns), dy(nt,ns), dz(nt,ns);

    for (int j=0; j<nt; j++) {
        for (int k=0; k<ns; k++) {
            dx(j,k)=(zetaf_r(j)-source_r(k));
            dy(j,k)=(zetaf_i(j)-source_i(k));
            dz(j,k)=dx(j,k)*dx(j,k)+dy(j,k)*dy(j,k)+delta(k)*delta(k);
            dx(j,k) = dx(j,k)/dz(j,k);
            dy(j,k) = dy(j,k)/dz(j,k);
        }
    }

    pot_r = -dx*distr/(2*ns*M_PI);
    pot_i = dy*distr/(2*ns*M_PI);


    /*
       gettimeofday(&time2, NULL); ttl2 = (time2.tv_sec * 1000000 + time2.tv_usec);
       diffl12 = (double)(ttl2 - ttl1)/1000000;
       printf("Elapsed %f seconds!\n", diffl12);
       */

    return ;
}
