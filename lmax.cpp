#include <cmath>
#include <cassert>
//#include "mpi.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "load-mols.h"
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

namespace linear {

//----------------------------------------------------------------------------//

static double* A(0);
static double* y(0);
static double* params(0);

//----------------------------------------------------------------------------//

void allocate()
{
//    A = new double[training_set.size()*model_type::nparams()];
//    y = new double[training_set.size()];

//    params = new double[model_type::nparams()];
}


} // namespace linear

////////////////////////////////////////////////////////////////////////////////

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

    // Get CVs

    for (size_t n = 0; n < training_set.size(); ++n) {
    
        training_set[n].xyz

    }

    
    
     
    
    
    
    
    



    std::cout << "Final parameters:" << std::endl;
//    for(size_t i = 0; i < num_nonlinear_params; ++i)
//        std::cout << all_atm_types[i]<< ' ' << s->x->data[i] << std::endl;
/*
    double err_L2(0), err_wL2(0), err_Linf(0);
    double err_L2_lo(0), err_Linf_lo(0), nlo(0);

    const char fm[] = "fit-output.txt";
    FILE *fs = fopen(fm, "w");
    fprintf(fs,"%6s %15s %15s %15s\n","index","model","training-set","delta");

    double emin = 1.0e+6;
    for (size_t n = 0; n < training_set.size(); ++n) {
        double mmm[model_type::nparams()];

///////////////////////////
        int nnlin = model_type::num_nonlinear_params;
        int ncoef = model_type::nparams();

        double tmp[nnlin+ncoef];

        for(size_t j=0; j<fit[0]->npar; ++j)
            tmp[j] = GA->parent[0][j];


        double data[nnlin];
        std::copy(tmp, tmp + nnlin , data);

        model.set_nonlinear_parameters(data);

        if (model.nonlinear_parameters_out_of_range())
          return 1.0e+6;

        double poly[ncoef];
        std::copy(tmp + nnlin, tmp + nnlin + ncoef, poly);

        model.set(poly);
        double chisq(0);

///////////////////////////

        double E_m = training_set[n].energy_total_com;
        double E_r = training_set[n].energy_total_ref;
//        double emin = 0.0;

//        const double E_model = model.value(training_set[n].xyz);
        const double delta = E_m - E_r;
        if (std::abs(delta) > err_Linf)
            err_Linf = std::abs(delta);

        err_L2 += delta*delta;
        err_wL2 += ts_weights[n]*delta*delta;

        if (E_r - E_min < E_range) {
            nlo += 1.0;
            err_L2_lo += delta*delta;
            if (std::abs(delta) > err_Linf_lo)
                err_Linf_lo = std::abs(delta);
        }

        fprintf(fs,"%06d %15.10f %15.10f %15.10f\n",n,E_m,E_r,delta);
    }

    fclose(fs);

    err_L2 /= training_set.size();
    err_wL2 /= training_set.size();

    err_L2_lo /= nlo;

    std::cout << "      err[L2] = " << std::sqrt(err_L2) << "\n"
              << "     err[wL2] = " << std::sqrt(err_wL2) << "\n"
              << "    err[Linf] = " << err_Linf << "\n"
              << "  err[L2,low] = " << err_L2_lo << "\n"
              << "err[Linf,low] = " << err_Linf_lo
              << std::endl;

    // output metadata
    std::cout << "\n"
              << "Degree polynomial = " << degree << "\n"
              << "Energy range      = " << E_range << "\n"
              << "Number of configs = " << training_set.size() << "\n"
              << std::endl;

*/

//    MPI_Finalize();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
