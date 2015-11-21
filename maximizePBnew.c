
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "./globals.txt"
#define Q (A+B)
#define DEBUG 1

/*************************************************************************/
/***********************  Function Declarations  *************************/
/*************************************************************************/

#include "./minroutinesPBnew.c"

/*************************************************************************/
/********************************  MAIN  *********************************/
/*************************************************************************/

int main(){

  int i, j, k, test, mode, mode1;
  double zA[1+Q][1+m], zB[1+Q][1+m], zavg[1+m];
  double movesize, gradsize, dlnLsize, maxsize, temp;
  double crit_move, crit_grad, crit_dlnL, x[1+m];
  double mlnlold, dmlnLold[1+m], H[1+m][1+m], K[1+m], S[1+m][1+m];
  double rc, a[1+m], aold[1+m], alpha0[1+m];
  double mlnl, dmlnL[1+m], alnLmax[1+m+1+1];
  FILE *convergence;

  getdata(zA, zB);

  for(j=0;j<=m;j=j+1){
    zavg[j] = 0.0;
  }
  for(i=1;i<=zA[0][0];i=i+1){
    for(j=1;j<=m;j=j+1){
      zavg[j] = zavg[j] + zA[i][j];
    }
  }
  for(j=1;j<=m;j=j+1){
    zavg[j] = zavg[j]/((double)zA[0][0]);
    alpha0[j] = 0.0;
  }

  crit_move = 0.001;
  crit_grad = 0.01;
  crit_dlnL = 0.0001;
  maxsize = 0.1;

  /*****  SCREEN FOR STARTING VALUES  *****/

  if(m == 1){
    printf("\n  GET INITIAL VALUES BY RANDOM NUMBERS");
    for(i=1;i<=16;i=i+1){
      a[1] = randomf(-2.0, 2.0);
      a[0] = -a[1]*zavg[j];
      mlnl = mlnL(a, zA, zB);
      if( mlnl < mlnL(alpha0, zA, zB) ){
        for(j=0;j<=1;j=j+1){
          alpha0[j] = a[j];
          printf("\n improvement %lf a: %.3f %.3f", 
                 mlnl, a[0], a[1]);
        }
      }
    }
  }else{
    printf("\n  GET INITIAL VALUES FROM BEST SUB-MODEL");
    convergence = fopen("initial.txt", "r");
    for(i=0;i<=m;i=i+1){
      fscanf(convergence, "%lf", &temp);
      alpha0[i] = temp;
    }
    fclose(convergence);
  }
  printf("\n  initial a: ");
  for(j=0;j<=m;j=j+1){
    a[j] = alpha0[j];
    printf(" %.3f ", a[j]);
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
    mlnl = grad(a, zA, zB, dmlnL);
    printf("\n dmlnL = %f %f", dmlnL[0], dmlnL[1]);
    if(i!=1){
      BFGS_update(x, dmlnL, dmlnLold, H);
    }
    diagonalize(H, S, K);
    dlnLsize = dmlnL - dmlnLold;
    gradsize = sqrt(dot(dmlnL, dmlnL));
    fprintf(convergence,"\n %d %f %f", i, gradsize, mlnl);
    printf("\n\n  gradients: %d  mlnL: %f  gradsize: %f", i, mlnl, gradsize);
    fflush(convergence);
    fflush(stdout);
    for(j=0;j<=m;j=j+1){
      x[j] = 0.0;
      for(k=0;k<=m;k=k+1){
	x[j] = x[j] - (dot(S[k],dmlnL)/K[k])*S[k][j];
      }
    }
    movesize = sqrt(dot(x,x));
    printf("\n  RAW MOVESIZE/(B): %f", movesize);
    if(movesize > maxsize){
      for(j=0;j<=m;j=j+1){
	x[j] = x[j]*maxsize/movesize;
      }
    }      
    printf("\n  NEW MOVESIZE/(B): %f  ALPHA: ", sqrt(dot(x,x)));
    for(j=0;j<=m;j=j+1){
      aold[j] = a[j];
      a[j] = a[j] + x[j];
      printf(" %.3f", a[j]);
    }
    printf("\n  CONVERGENCE TESTING");
    printf("\n   movesize:tolerance  %f::%f", movesize, crit_move);
    printf("\n   gradsize:tolerance  %f::%f", gradsize, crit_grad);
    printf("\n   deltaE:|tolerance|  %f::%f", dlnLsize, crit_dlnL);
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
      printf("\n  CONVERGED ");
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

  return 0;
}


