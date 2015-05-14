/* This is the computational routine of mex file main routine ByFreeAccurate_singleSheet.cpp
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
void ByFreeAccurate_singleSheet_routine(size_t nt, double *ztreal, double *ztimag, size_t ns, double *zetaf_r, double *zetaf_i, double *potreal, double *potimag, double *L, double *gammaf)
{
    timeval time1, time2;
    double tt1, tt2, diff12;
    gettimeofday(&time1, NULL); tt1 = (time1.tv_sec * 1000000 + time1.tv_usec);

    VectorXd zztreal(nt), zztimag(nt);
    VectorXd ddistr(ns);
    VectorXd zzetaf_r(ns), zzetaf_i(ns);

    /* source is all the free sheets contaminated */
    /* read free sheets from matlab to c++ */
    for (int i = 0; i < nt; ++i){
    zztreal(i) = ztreal[i];
    zztimag(i) = ztimag[i];
    }
    vector<Wing> W(10);
    W[0].zetaf_r = zztreal.block<10,1>(0,0);
    W[0].zetaf_i = zztimag.block<10,1>(0,0);

    W[1].zetaf_r = zztreal.block<10,1>(10,0);
    W[1].zetaf_i = zztimag.block<10,1>(10,0);

    for (int i = 0; i < ns; ++i){
    zzetaf_r(i) = zetaf_r[i];
    zzetaf_i(i) = zetaf_i[i];
    }
    for (int i = 0; i < ns; ++i){
        ddistr(i) = gammaf[i];
    }

    /* Computing */
    VectorXd pot_r, pot_i;

    /*
*/
    omp_set_num_threads(8);
#pragma omp parallel for
    for (int it = 0; it < 8; ++it){
/*    printf("Thread %d is starting...\n",omp_get_thread_num());
 */
    byfree(zztreal, zztimag, zzetaf_r, zzetaf_i, ddistr, pot_r, pot_i);
    }

    /* Convert output from eigen vector to array */
    for (int j=0; j<nt; j++) {
        potreal[j] = pot_r(j);
        potimag[j] = pot_i(j);
    }

    gettimeofday(&time2, NULL); tt2 = (time2.tv_sec * 1000000 + time2.tv_usec);
    diff12 = (double)(tt2 - tt1)/1000000;
    printf("Elapsed %f seconds!\n", diff12);
}





/* Implementation of byfree */
void byfree(VectorXd target_r, VectorXd target_i, VectorXd zetaf_r, VectorXd zetaf_i, VectorXd gammaf, VectorXd& pot_r, VectorXd& pot_i){
    timeval time1, time2;
    double tt1, tt2, diff12;
    gettimeofday(&time1, NULL); tt1 = (time1.tv_sec * 1000000 + time1.tv_usec);

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

    gettimeofday(&time2, NULL); tt2 = (time2.tv_sec * 1000000 + time2.tv_usec);
    diff12 = (double)(tt2 - tt1)/1000000;
    printf("Elapsed %f seconds!\n", diff12);

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

  */

