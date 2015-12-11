#ifndef MAXIMIZE_H
#define MAXIMIZE_H

#include <cstdlib>
#include <vector>
//#include "dsyev.h"

double *x;//        = new double [1+m];
double *zavg;//     = new double [1+m];
double *dmlnLold;// = new double [1+m];
double *K;//        = new double [1+m];
double *a;//        = new double [1+m];
double *aold;//     = new double [1+m];
double *alpha0;//   = new double [1+m];
double *dmlnL;//    = new double [1+m];
double *alnLmax;//  = new double [1+m+1+1];
    
double **zA;//    = new double* [1+Q];
double **zB;//    = new double* [1+Q];
  
double **H;//     = new double* [1+m];
double **S;//     = new double* [1+m];

int maximize(int NA, int NB, int m);

#endif // MAXIMIZE_H
