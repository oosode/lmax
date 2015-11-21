#include <cmath>
#include <cassert>
//#include "mpi.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

#include "load-mols.h"
#include "load-data.h"
#include "constants.h"

//#define VERBOSE 1

#define DEBUG 1 
#define MAX_LENGTH   64

#ifdef DEBUG
#define PR(x) std::cout << #x << ": " << (x) << std::endl;
#else
#define PR(x)
#endif /* DEBUG */

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//
const int degree = 3;
const double tol = 1e-3;

//typedef offdiag::offdiag_v1<degree> model_type;
//static model_type model;
    
#ifdef RIDGE_REGRESSION
const double alpha = 0.0005;
#endif
const double E_range = 25/constants::Eh_kcalmol; // kcal/mol

static std::vector<shootpt::mol> training_set;
static double*                   ts_weights = 0;

    
////////////////////////////////////////////////////////////////////////////////

double cv1(const int n)
{

    shootpt::mol ts = training_set[n];
    double cv = 0;
    double r0 = 2.38;

    int beta[3] = {5867,5868,5869}; // 5868, 5869, 5870
    int gammap  = 5870; // 5871

    for (size_t i=0; i<3; ++i) {

        double d2 = 0;
        double  d = 0;

        for (size_t l=0; l<3; ++l) {

            d = ts.xyz[gammap*3+l]-ts.xyz[beta[i]*3+l];
	    d2 += d*d;
        
        }

	d = sqrt (d2);

        double num = 1 - pow(d/r0,6.0);
        double den = 1 - pow(d/r0,12.0);

        cv += num/den;        

    }
    
    return cv;
}

double cv2(const int n)
{

    shootpt::mol ts = training_set[n];
    double cv = 0;
    double r0 = 2.38;

    int waterstart = 5875; //5876
    int waterend   = 55749; // 55750

    
    int gammao[3] = {5871,5872,5873}; // 5872, 5873, 5874
    int gammap  = 5870; // 5871

    for (size_t i=0; i<3; ++i) {

        double d2 = 0;
        double  d = 0;

        for (size_t l=0; l<3; ++l) {

            d = ts.xyz[gammap*3+l]-ts.xyz[gammao[i]*3+l];
            d2 += d*d;

	}
	
	d = sqrt (d2);

	double num = 1 - pow(d/r0,6.0);
	double den = 1 - pow(d/r0,12.0);

	cv += num/den;

    }

    for (size_t i=waterstart; i<waterend; i+=3) {

        double d2 = 0;
        double  d = 0;

        for (size_t l=0; l<3; ++l) {

            d = ts.xyz[gammap*3+l]-ts.xyz[i*3+l];
            d2 += d*d;

        }

        if (d2>64) continue;

        d = sqrt (d2);

        double num = 1 - pow(d/r0,6.0);
        double den = 1 - pow(d/r0,12.0);

        cv += num/den;

    }

    return cv;
}


} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

//    MPI_Init(&argc, &argv);

//    MPI_Comm_size(MPI_COMM_WORLD, &(global.nproc));
//    MPI_Comm_rank(MPI_COMM_WORLD, &(global.iproc));

    //
    // initialize the model and load training sets
    //

    if (argc < 2) {
        std::cerr << "usage: lmax.x x ..."
                  << std::endl;
        return 0;
    }

//    std::cout << "\n<><><> model type = '" << model_type::name() << " <><><>'\n";

    ++argv;

    try {
        while (--argc != 0) {
            size_t nd = shootpt::load_mols(training_set);
            std::cout << "'" << *(argv++) << "' : "
                      << nd << " mols" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    //Get accepted/rejected 
    char basin_file[64];
    char acc_file[64];

    sprintf(acc_file, "acc.txt");
    sprintf(basin_file, "basin_evals.txt");

    shootpt::load_acc(training_set, acc_file);
    shootpt::load_basin(training_set, basin_file);

    // Get CVs
    int M = 2; // number of CVs
    double qA[M][training_set.size()], qB[M][training_set.size()];
    double zA[M][training_set.size()], zB[M][training_set.size()];
    double qmax[M], qmin[M], qspan[M];
    double qmaxA[M], qmaxB[M], qminA[M], qminB[M];

    size_t NA = 0;
    size_t NB = 0;

    for (size_t n = 0; n < training_set.size(); ++n) {
 
	shootpt::mol ts = training_set[n];

        if (ts.conclusive == 1) {

	    if (ts.haf == 1 && ts.hbb == 1) {
		qA[0][NA] = cv1(n);
		std::cout << "qA " << n << " " << qA[0][NA] << std::endl;
		qA[1][NA] = cv2(n);
		NA++;
	    } else if (ts.hbf == 1 && ts.hab == 1) {
		qB[0][NB] = cv1(n);
	        std::cout << "qB " << n << " " << qB[0][NB] << std::endl;
		qB[1][NB] = cv2(n);
		NB++;
	    }
	}
    }

 
    for (size_t i = 0; i < M; ++i) {
//        if (ts.haf == 1 && ts.hbb == 1) {
            qminA[i] = qA[i][0];
	    qmaxA[i] = qA[i][0];
	    for (size_t j = 0; j < NA; ++j) {
                if (qA[i][j] < qminA[i]) {
                    qminA[i] = qA[i][j];
                }
                if (qA[i][j] > qmaxA[i]) {
                    qmaxA[i] = qA[i][j];
		}
            }
//	} else if (ts.hbf == 1 && ts.hab == 1) {
            qminB[i] = qB[i][0];
	    qmaxB[i] = qB[i][0];
            for (size_t j = 0; j < NB; ++j) {
                if (qB[i][j] < qminB[i]) {
                    qminB[i] = qB[i][j];
                }
                if (qB[i][j] > qmaxB[i]) {
                    qmaxB[i] = qB[i][j];
                }
            }
//	}
    }

    printf("\nCV qmin  qmax  qspan\n");
    for (size_t i = 0; i < M; ++i) {
        if (qmaxA[i] > qmaxB[i]) {
            qmax[i] = qmaxA[i];
        } else {
            qmax[i] = qmaxB[i];
        }
        if (qminA[i] < qminB[i]) {
            qmin[i] = qminA[i];
        } else{
            qmin[i] = qminB[i];
        }
        qspan[i] = qmax[i] - qmin[i];
        printf("\n %zu %.3f %.3f %.3f ", i, qmin[i], qmax[i], qspan[i]);
    }

    printf("\n\n  REDUCED VARIABLES: Z = (Q-Qmin)/(Qmax-Qmin) \n");
    printf("\n        Z IN [0, 1]: Q = Z(Qmax-Qmin)+Qmin\n");    
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < NA; ++j) {
            zA[i][j] = (qA[i][j]-qmin[i])/(qmax[i]-qmin[i]);
        }
        for (size_t j = 0; j < NB; ++j) {
            zB[i][j] = (qB[i][j]-qmin[i])/(qmax[i]-qmin[i]);
        }
    }  
    
    FILE *fa = fopen("A.txt", "w");

    fprintf(fa,"%d\n",NA);
    for(int j=0; j<NA; j++){
      fprintf(fa, "\n");
      for(int i=0; i<M; i++){
        fprintf(fa, "%f ", zA[i][j]);
      }
    }
    fclose(fa);

    FILE *fb = fopen("B.txt", "w");

    fprintf(fb, "%d\n", NB);
    for(int j=0; j<NB; j++){
      fprintf(fb, "\n");
      for(int i=0; i<M; i++){
        fprintf(fb, "%f ", zB[i][j]);
      }
    }
    fclose(fb);

    printf("\n\n  BIC = LN(N)/2 = %.3f\n", 0.5*log((double)(NA+NB)));    

    printf("\n\n PREPARE AND COMPILE MAXIMIZATION ROUTINE \n");

    FILE *fm = fopen("globals.txt","w");
    fprintf(fm, "\n#define A %d", NA);
    fprintf(fm, "\n#define B %d", NB);
    fprintf(fm, "\n#define m %d\n", 1);
    fclose(fm);

    system("gcc -lm maximizePBnew.c");
    system("cp a.out max.exe");
    
    int imax;
    double rc, alpha[M], lnL, temp;
    double best1[M][1+1+1], best2[M][M][1+2+1];
    double best3[M][M][M][1+3+1];

    int m = 1;

    for(int i=0; i<m; i++){
      FILE *zlist = fopen("zlist.txt", "w");
      fprintf(zlist, "\n %d %d\n", NA, NB);
      for(int j=0; j<NA; j++){
        fprintf(zlist, "\n %f ", zA[i][j]);
      }
      for(int j=0; j<NB; j++){
        fprintf(zlist, "\n %f ", zB[i][j]);
      }
      fclose(zlist);
      system("./max.exe > a1.txt");
      fflush(stdout);
      FILE *answer = fopen("answer.txt", "r");
      fscanf(answer, "%lf ", &temp);
      lnL = temp;
      best1[i][m] = lnL;
      if (i==1) {
        system("cp zlist.txt z-1.txt");
        imax = 1;
      }else{
        if(best1[i][m] > best1[imax][m]){
  	  system("cp zlist.txt z-1.txt");
	  imax = i;
        }
      }
      printf("\n rxncoor %d %f", i, lnL);
      for (int j=0; j<m; j++){
        fscanf(answer, "%lf ", &temp);
        alpha[j] = temp;
        printf(" %.3f", alpha[j]);
        best1[i][j] = alpha[j];
      }
      fclose(answer);
    }
    printf("\n ");
    printf("\n best rxncoor %d %f", imax, best1[imax][m]);
    for(int j=0; j<m; j++){
      printf(" %.3f", best1[imax][j]);
    }


    std::cout << "Final parameters:" << std::endl;


    return 0;
}

////////////////////////////////////////////////////////////////////////////////
