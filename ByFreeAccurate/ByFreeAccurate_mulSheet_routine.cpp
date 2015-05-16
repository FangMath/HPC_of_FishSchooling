/* This is the computational routine of mex file main routine ByFreeAccurate_mulSheet.cpp
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

void byfree(VectorXd target_r, VectorXd target_i, VectorXd zetaf_r, VectorXd zetaf_i, VectorXd gammaf, VectorXd& pot_r, VectorXd& pot_i);

/* The computational routine */
void ByFreeAccurate_mulSheet_routine(size_t nt, double *ztreal, double *ztimag, size_t ns, double *zetaf_r, double *zetaf_i, double *potreal, double *potimag, size_t Nw, double *L, double *gammaf)
{
    /*
    cout << "Nw=" << Nw << endl;
    timeval time1, time2;
    double tt1, tt2, diff12;
    gettimeofday(&time1, NULL); tt1 = (time1.tv_sec * 1000000 + time1.tv_usec);
    */

    VectorXd zztreal(nt), zztimag(nt);
    VectorXd ggammaf(ns);
    VectorXd zzetaf_r(ns), zzetaf_i(ns);

    /* read data from array to eigen vectors */
    for (int i = 0; i < nt; ++i){
    zztreal(i) = ztreal[i];
    zztimag(i) = ztimag[i];
    }

    for (int i = 0; i < ns; ++i){
    zzetaf_r(i) = zetaf_r[i];
    zzetaf_i(i) = zetaf_i[i];
    }
    for (int i = 0; i < ns; ++i){
        ggammaf(i) = gammaf[i];
    }

    /* read free sheets from array to c++ class */
    vector<Wing> W(Nw);
    W[0].zetaf_r = zzetaf_r.head(L[0]);
    W[0].zetaf_i = zzetaf_i.head(L[0]);
    W[0].gammaf = ggammaf.head(L[0]);

    for (int ib=1; ib<Nw; ++ib)
    {
    W[ib].zetaf_r = zzetaf_r.segment(L[ib-1], L[ib]-L[ib-1]);
    W[ib].zetaf_i = zzetaf_i.segment(L[ib-1], L[ib]-L[ib-1]);
    W[ib].gammaf = ggammaf.segment(L[ib-1], L[ib]-L[ib-1]);
    }


    /* Computing */
    VectorXd pot_r = VectorXd::Zero(nt,1), pot_i = VectorXd::Zero(nt,1);

cout << "max threads=" << omp_get_max_threads() << endl;
int Np = Nw;
omp_set_num_threads(Np);
#pragma omp parallel
    {
        VectorXd ppot_r = VectorXd::Zero(nt,1), ppot_i = VectorXd::Zero(nt,1);

#pragma omp for nowait 
        for (int ib = 0; ib < Nw; ++ib){
        VectorXd tpot_r = VectorXd::Zero(nt,1), tpot_i = VectorXd::Zero(nt,1);
            /*    
             * printf("Thread %d is starting...\n",omp_get_thread_num());
            */
            byfree(zztreal, zztimag, W[ib].zetaf_r, W[ib].zetaf_i, W[ib].gammaf, tpot_r, tpot_i);
            ppot_r += tpot_r;
            ppot_i += tpot_i;
        }

#pragma omp critical
        {
            pot_r += ppot_r;
            pot_i += ppot_i;
        }

    }


    /* Convert output from eigen vector to array */
    for (int j=0; j<nt; j++) {
        potreal[j] = pot_r(j);
        potimag[j] = pot_i(j);
    }

    /*
    gettimeofday(&time2, NULL); tt2 = (time2.tv_sec * 1000000 + time2.tv_usec);
    diff12 = (double)(tt2 - tt1)/1000000;
    printf("Elapsed %f seconds!\n", diff12);
    */
}





/* Implementation of byfree */
void byfree(VectorXd target_r, VectorXd target_i, VectorXd zetaf_r, VectorXd zetaf_i, VectorXd gammaf, VectorXd& pot_r, VectorXd& pot_i){

    /*
    timeval time1, time2;
    double ttl1, ttl2, diffl12;
    gettimeofday(&time1, NULL); ttl1 = (time1.tv_sec * 1000000 + time1.tv_usec);
    */

    int nt = target_r.size(), ns = zetaf_i.size();

    MatrixXd dx(nt,ns), dy(nt,ns), dz(nt,ns), F_r(nt,ns-1), F_i(nt,ns-1);
    VectorXd distr, distr_r(ns-1), distr_i(ns-1);

    distr = -(gammaf.tail(ns-1) - gammaf.head(ns-1));
        for (int k=0; k<ns-1; k++) {
            double dxnorm = (zetaf_r(k+1) - zetaf_r(k))*(zetaf_r(k+1) - zetaf_r(k)) + (zetaf_i(k+1) - zetaf_i(k))*(zetaf_i(k+1) - zetaf_i(k));
            distr_r(k) = (zetaf_r(k+1) - zetaf_r(k))*distr(k)/dxnorm;
            distr_i(k) = -(zetaf_i(k+1) - zetaf_i(k))*distr(k)/dxnorm;

        }


    for (int j=0; j<nt; j++) {
        for (int k=0; k<ns; k++) {
            dx(j,k)=(target_r(j)-zetaf_r(k));
            dy(j,k)=(target_i(j)-zetaf_i(k));
            dz(j,k)=dx(j,k)*dx(j,k)+dy(j,k)*dy(j,k);
            }
        }

    for (int j=0; j<nt; j++) {
        for (int k=0; k<ns-1; k++) {
            F_r(j,k) = log(dz(j,k)/dz(j,k+1))/2;
            F_i(j,k) = atan2((-dx(j,k)*dy(j,k+1)+dy(j,k)*dx(j,k+1)),(dx(j,k)*dx(j,k+1)+dy(j,k)*dy(j,k+1)));
        }
    }
    
    pot_r = F_r*distr_r - F_i*distr_i;
    pot_i = F_i*distr_r + F_r*distr_i;;

    /*
    gettimeofday(&time2, NULL); ttl2 = (time2.tv_sec * 1000000 + time2.tv_usec);
    diffl12 = (double)(ttl2 - ttl1)/1000000;
    printf("Elapsed %f seconds!\n", diffl12);
    */

    return ;
}


  /*
#pragma omp parallel for
    W[0].zetaf_r = ddistr;

    cout << W[1].K << endl;
    cout << "ff" << endl;
    #pragma omp parallel 
    {
       int  pn = omp_get_num_threads();
    if (omp_get_thread_num() == 0) {
         printf("Number of threads = %d\n", omp_get_num_threads()); }
  printf("Thread %d is starting...\n",omp_get_thread_num());
#pragma omp parallel for
    } 
    printf("Max of threads = %d\n", omp_get_max_threads()); 
    omp_set_num_threads(8);
    #pragma omp parallel 
    {
       int  pn = omp_get_num_threads();
    if (omp_get_thread_num() == 0) {
         printf("Number of threads = %d\n", omp_get_num_threads()); }
  printf("Thread %d is starting...\n",omp_get_thread_num());
    } 

#pragma omp parallel
    {
        VectorXd ppot_r, ppot_i;

#pragma omp for nowait 
        for (int ib = 0; ib < Nw; ++ib){
            byfree(zztreal, zztimag, W[ib].zetaf_r, W[ib].zetaf_i, W[ib].gammaf, ppot_r, ppot_i);
        }

#pragma omp critical
        {
            pot_r += ppot_r;
            pot_i += ppot_i;
        }

    }

  */

