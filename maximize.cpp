#include <cassert>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//#include "globals.txt"
//#define Q (A+B)
#define DEBUG 0

/*************************************************************************/
/***********************  Function Declarations  *************************/
/*************************************************************************/

#include "minroutines.h"

/*************************************************************************/
/********************************  MAIN  *********************************/
/*************************************************************************/

int maximize(int NA, int NB, int m) {

  int Q = NA+NB;
  
  int i, j, k, test, mode, mode1;
  double rc, mlnl;
//  double zA[1+Q][1+m], zB[1+Q][1+m], 
  double movesize, gradsize, dlnLsize, maxsize, temp;
  double crit_move, crit_grad, crit_dlnL;
  double mlnlold;

  double *x        = new double [1+m];
  double *zavg     = new double [1+m];
  double *dmlnLold = new double [1+m];
  double *K        = new double [1+m];
  double *a        = new double [1+m];
  double *aold     = new double [1+m];
  double *alpha0   = new double [1+m];
  double *dmlnL    = new double [1+m];
  double *alnLmax  = new double [1+m+1+1];

//  double **zA;
//  double **zB;
  double **zA    = new double*[1+Q];
  for (int l=0; l<=Q; ++l) {
    zA[l] = new double[1+m];
  }
  double **zB    = new double*[1+Q];
  for (int l=0; l<=Q; ++l) {
    zB[l] = new double[1+m];
  }
  zA[0][0]=NA;
  zB[0][0]=NB;


  double **H     = new double*[1+m];
  double **S     = new double*[1+m];
  for (i=0; i<=m; ++i) {
    H[i] = new double[1+m];
    S[i] = new double[1+m];
  }

  FILE *convergence;

//  printf("%d %d %f %f %d -- values in maximize\n",NA,NB,zA[0][0],zB[0][0],m);
  getdata(zA,zB,NA,NB,m);

  for(j=0;j<=m;j=j+1){
    zavg[j] = 0.0;
  }
  for(i=1;i<=(int)zA[0][0];i=i+1){
    for(j=1;j<=m;j=j+1){
      zavg[j] = zavg[j] + zA[i][j];
    }
  }
  for(j=1;j<=m;j=j+1){
    zavg[j] = zavg[j]/((double)zA[0][0]);
    alpha0[j] = 0.0;
  }
  //exit(0);

  crit_move = 0.001;
  crit_grad = 0.01;
  crit_dlnL = 0.0001;
  maxsize   = 0.1;

  /*****  SCREEN FOR STARTING VALUES  *****/
  j=1;
  if(m == 1) {
    if (DEBUG==1) printf("\n  GET INITIAL VALUES BY RANDOM NUMBERS");
    for(i=1;i<=16;i=i+1){
      a[1] = randomf(-2.0, 2.0);
      a[0] = -a[1]*zavg[1];
//      a[0] = -a[1]*zavg[j];
//      printf("%f %f %f %d -- nnnn\n",a[0],a[1],zavg[1],j);
      mlnl = mlnL(a, zA, zB, m);
      if( mlnl < mlnL(alpha0, zA, zB, m) ){
        for(j=0;j<=1;j=j+1){
          alpha0[j] = a[j];
          if (DEBUG==1) printf("\n improvement %lf a: %.3f %.3f", 
                 mlnl, a[0], a[1]);
        }
      }
    }
  } else {
    if (DEBUG==1) printf("\n  GET INITIAL VALUES FROM BEST SUB-MODEL");
    convergence = fopen("initial.txt", "r");
    for(i=0;i<=m;i=i+1){
      fscanf(convergence, "%lf", &temp);
      alpha0[i] = temp;
    }
    fclose(convergence);
  }
//  exit(0);
  if (DEBUG==1) printf("\n  initial a: ");
  for(j=0;j<=m;j=j+1){
    a[j] = alpha0[j];
    if (DEBUG==1) printf(" %.3f ", a[j]);
    for(k=0;k<=m;k=k+1){
      H[j][k] = 0.0;
    }
    H[j][j] = 1.0 + randomf(0.0, 0.1);
  }
  
  convergence = fopen("convergence.txt", "w");
  test = 0;
  for(i=1;test==0;i=i+1){
    if(i != 1){
      mlnlold = mlnl;
      for(j=0;j<=m;j=j+1){
	dmlnLold[j] = dmlnL[j];
      }
    }
    mlnl = grad(a, zA, zB, dmlnL, m);
    if (DEBUG==1) printf("\n dmlnL = %f %f", dmlnL[0], dmlnL[1]);
    if(i!=1){
      BFGS_update(x, dmlnL, dmlnLold, H, m);
    }
    diagonalize(H, S, K, m);
    dlnLsize = dmlnL - dmlnLold;
    gradsize = sqrt(dot(dmlnL, dmlnL,m));
    fprintf(convergence,"\n %d %f %f", i, gradsize, mlnl);
    if (DEBUG==1) printf("\n\n  gradients: %d  mlnL: %f  gradsize: %f", i, mlnl, gradsize);
    fflush(convergence);
    fflush(stdout);
    for(j=0;j<=m;j=j+1){
      x[j] = 0.0;
      for(k=0;k<=m;k=k+1){
	x[j] = x[j] - (dot(S[k],dmlnL,m)/K[k])*S[k][j];
      }
    }
    movesize = sqrt(dot(x,x,m));
    if (DEBUG==1) printf("\n  RAW MOVESIZE/(B): %f", movesize);
    if(movesize > maxsize){
      for(j=0;j<=m;j=j+1){
	x[j] = x[j]*maxsize/movesize;
      }
    }      
    if (DEBUG==1) printf("\n  NEW MOVESIZE/(B): %f  ALPHA: ", sqrt(dot(x,x,m)));
    for(j=0;j<=m;j=j+1){
      aold[j] = a[j];
      a[j] = a[j] + x[j];
      if (DEBUG==1) printf(" %.3f", a[j]);
    }
   
    if (DEBUG==1) {
      printf("\n  CONVERGENCE TESTING");
      printf("\n   movesize:tolerance  %f::%f", movesize, crit_move);
      printf("\n   gradsize:tolerance  %f::%f", gradsize, crit_grad);
      printf("\n   deltaE:|tolerance|  %f::%f", dlnLsize, crit_dlnL);
    }
    if(movesize < crit_move && gradsize < crit_grad){
      test = 1;
    }
    if(movesize < crit_move && fabs(dlnLsize) < crit_dlnL){
      test = 1;
    }
    if(fabs(dlnLsize) < crit_dlnL && gradsize < crit_grad){
      test = 1;
    }
    if(test == 1){
      if (DEBUG==1) printf("\n  CONVERGED ");
      for(j=0;j<=m;j=j+1){
	alnLmax[j] = a[j];
      }
      alnLmax[m+1] = -mlnl;
    }
  }
  fclose(convergence);

  convergence = fopen("answer.txt", "w");
  fprintf(convergence, "%f ", alnLmax[m+1]);
  for(j=0;j<=m;j=j+1){
    fprintf(convergence, "%f ", alnLmax[j]);
  }
  fclose(convergence);
  printf("\n");
//  exit(0);
  return 0;
}


