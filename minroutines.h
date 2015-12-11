#ifndef MINROUTINES_H
#define MINROUTINES_H

//double rxn_coor(double *a, double *z);
//double mlnL(double *a, double **zA, double **zB);
//double dot(double [1+m], double [1+m]);
//double normalize(double [1+m]);
//void clone(double [1+m], double [1+m]);
//double sech(double);
//double pA(double);
//void getdata(double [1+Q][1+m], double [1+Q][1+m]);
//void diagonalize(double [1+m][1+m], double [1+m][1+m], double [1+m]);
//double randomf(double, double);
//void projectfrommatrix(double [1+m], double [1+m][1+m]);
//void BFGS_update(double [1+m], double [1+m], 
//		 double [1+m], double [1+m][1+m]);
//double grad(double *a, double **zA, double **zB, double *dmlnL);



//double grad(double a[1+m], double zA[1+Q][1+m], 
//	    double zB[1+Q][1+m], double dmlnL[1+m]){


double dot(double *v, double *u, int m) {
    int k;
    double dp;
    dp = 0.0;
    for(k=0;k<=m;k=k+1){
        dp = dp+v[k]*u[k];
    }
    return dp;
}

double rxn_coor(double *a, double *z, int m) {
    double rc;
    int i;
    rc = a[0];
    for(i=1;i<=m;i=i+1){
        rc = rc + z[i]*a[i];
    }
    return rc;
}

double mlnL(double *a, double **zA, double **zB, int m) {
    int i, j;
    double lnl, rk;
//    printf("%f lnl",(double)(zA[0][0]+zB[0][0]));
    lnl = ((double)(zA[0][0]+zB[0][0]))*log(0.5);
    for(i=1;i<=zA[0][0];i=i+1){
        rk = rxn_coor(a, zA[i],m);
        lnl = lnl + log(1.0 - tanh(rk));
    }
    for(i=1;i<=zB[0][0];i=i+1){
        rk = rxn_coor(a, zB[i],m);
        lnl = lnl + log(1.0 + tanh(rk));
    }
    return (-lnl);
}

double normalize(double *x, int m) {
    
    double numer;
    int j;
    numer = 0.0;
    for(j=0;j<=m;j=j+1){
        numer = numer + x[j]*x[j];
    }
    numer = sqrt(numer);
    for(j=0;j<=m;j=j+1){
        x[j] = x[j]/numer;
    }
    return numer;
}

double grad(double *a, double **zA, double **zB, double *dmlnL, int m) {

  int i, j, k;
  double lnl, rk, C;

  lnl = -mlnL(a, zA, zB, m);
  for(i=0;i<=m;i=i+1){
    dmlnL[i] = 0.0;
  }
  for(i=1;i<=zA[0][0];i=i+1){
    rk = rxn_coor(a, zA[i],m);
    C = 1.0+tanh(rk);
    dmlnL[0] = dmlnL[0] + C;
    for(j=1;j<=m;j=j+1){
      dmlnL[j] = dmlnL[j] + C*zA[i][j];
    }
  }
  for(i=1;i<=zB[0][0];i=i+1){
    rk = rxn_coor(a, zB[i],m);
    C = 1.0-tanh(rk);
    dmlnL[0] = dmlnL[0] - C;
    for(j=1;j<=m;j=j+1){
      dmlnL[j] = dmlnL[j] - C*zB[i][j];
    }
  }
  return (-lnl);
}

double randomf(double a, double b) {
    
    double aa;
    aa = ((double)( rand()%10001))/10000.0;
    return a+(b-a)*aa;
    
}

void projectfrommatrix(double *vector, double **hmwc, int m) {

    int i,j,k;
    double p[1+m][1+m];
    double K[1+m][1+m];
    
    for(i=0;i<=m;i=i+1){
        for(j=0;j<=m;j=j+1){
            if(i == j){
                p[i][j] = 1.0 - vector[i]*vector[j];
            }
            else{
                p[i][j] = - vector[i]*vector[j];
            }
        }
    }
    for(i=0;i<=m;i=i+1){
        for(j=0;j<=m;j=j+1){
            K[i][j] = 0.0;
            for(k=0;k<=m;k=k+1){
                K[i][j] = K[i][j] + hmwc[i][k]*p[k][j];
            }
        }
    }
    for(i=0;i<=m;i=i+1){
        for(j=0;j<=m;j=j+1){
            hmwc[i][j] = 0.0;
            for(k=0;k<=m;k=k+1){
                hmwc[i][j] = hmwc[i][j] + p[i][k]*K[k][j];
            }
        }
    }
}

void diagonalize(double **hmwc, double **smwc, double *w2, int m) {
  
  int i, j, k, ndiag, test;
  double diff, tolerance, relerror, evalue;
//  double K[1+m][1+m], x[3][1+m];
  double x[3][1+m];

  double **K = new double *[1+m];
  for (i=0; i<=m; ++i) {
    K[i] = new double [1+m];
  }

  ndiag = 1000;

  for(j=0;j<=m;j=j+1){
    for(k=0;k<=m;k=k+1){
      K[j][k] = hmwc[j][k];
    }
  }
  tolerance = .000001;
  if(DEBUG == 1){
    printf("\n  DIAGONALIZING:");
    fflush(stdout);
  }
  for(i=0;i<=m;i=i+1){
    for(j=0;j<=m;j=j+1){
      x[0][j]=randomf(-1.0,1.0);
      x[1][j] = 0.0;
      x[2][j] = 0.0;
    }
    normalize(x[0],m);
    test = 0;
    for(k=1;test == 0;k=k+1){
      for(j=0;j<=m;j=j+1){
        x[k%3][j]=dot(K[j],x[(k-1)%3],m);
      }
      normalize(x[k%3],m);
      relerror = 0.0;
      for(j=0;j<=m;j=j+1){
        diff = x[k%3][j]-x[(k-2)%3][j];
        relerror = relerror + fabs(diff);
      }
      if((relerror < tolerance && k > 100) || (k == ndiag)){
        if(DEBUG == 1){
          printf(".");
          if(k == ndiag){
            printf("\n  itmax reached:");
            printf("    relative error in evect %d: %f",i,relerror);
          }
        }
        for(j=0;j<=m;j=j+1){
          x[(k+1)%3][j] = dot(K[j],x[k%3],m);
        }
        w2[i] = dot(x[k%3],x[(k+1)%3],m);
        if(DEBUG == 1){
          printf(" (%.2f) ",w2[i]);
        }
        for(j=0;j<=m;j=j+1){
          smwc[i][j] = x[k%3][j];
        }
        projectfrommatrix(smwc[i],K,m);
        test = 1;
      }
    }
  }
// BUGBUGBUGBUG
//  exit(0);
}

void BFGS_update(double *dx, double *g_new, double *g_old, double **h, int m){

  int i, j, k;
  double Hdx[1+m], y[1+m], dxy;
  double r[1+m], theta, dxHdx, dxr;
  
  if (DEBUG==1) {
    printf("\n  DAMPED BFGS HESSIAN UPDATE p.541 NOCEDAL & WRIGHT");
    printf("\n    SECANT");
  }
  for(j=0;j<=m;j=j+1){
    y[j] = g_new[j] - g_old[j];
  }
  
  if (DEBUG==1) printf("  H.dx:");
  for(i=0;i<=m;i=i+1){
    Hdx[i] = dot(h[i], dx, m);
    if (DEBUG==1) printf(" %.2f",Hdx[i]);
  }
  dxHdx = dot(dx, Hdx, m);
  dxy = dot(dx, y, m);
  if( dxy > 0.2*dxHdx ){
    theta = 1.0;
  }else{
    theta = 0.8*dxHdx/(dxHdx-dxy);
  }
  if (DEBUG==1) printf("  theta: %.3f", theta);
  for(i=0;i<=m;i=i+1){
    r[i] = theta*y[i]+(1.0-theta)*Hdx[i];
  }
  dxr = dot(r, dx, m);
  for(i=0;i<=m;i=i+1){
    for(j=0;j<=m;j=j+1){
      h[i][j] = h[i][j] - Hdx[i]*Hdx[j]/dxHdx + r[i]*r[j]/dxr;
    }
  }
}

void getdata(double **zA, double **zB, int NA, int NB, int m) {

  FILE *zlist;
  int i, j, k, itemp;
  double temp;

  zlist = fopen("zlist.txt", "r");
  fscanf(zlist, "%d ", &itemp);
  if(itemp != NA){
    printf("\n  DISCREPANCY IN EXPECTED & WRITTEN LENGTHS OF A-LIST");
  }
  fscanf(zlist, "%d ", &itemp);
  if(itemp != NB){
    printf("\n  DISCREPANCY IN EXPECTED & WRITTEN LENGTHS OF B-LIST");
  }
  for(j=1;j<=(int)zA[0][0];j=j+1){
    for(k=1;k<=m;k=k+1){
      fscanf(zlist, "%lf ", &temp);
      zA[j][k] = temp;
    }
  } 
  for(j=1;j<=(int)zB[0][0];j=j+1){
    for(k=1;k<=m;k=k+1){
      fscanf(zlist, "%lf ", &temp);
      zB[j][k] = temp;
    }
  }
}

double sech(double x) {

  double temp;
  temp = cosh(x);
  return (1.0/temp);

}

void clone(double *x, double *y, int m) {

  int i;
  for(i=0;i<=m;i=i+1){
    y[i] = x[i];
  }

}


#endif // MINROUTINES_H
