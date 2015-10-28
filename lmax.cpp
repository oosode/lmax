#include <cmath>
#include <cassert>
#include "mpi.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <gsl/gsl_multimin.h>

#define RIDGE_REGRESSION 1

#ifdef RIDGE_REGRESSION
#   include "rwlsq.h"
#else
#   include "wlsq.h"
#endif

#include "offdiag-v1.h"
#include "x2b-h3o-h2o-v1.h"
#include "x1b-hydronium-v1.h"
#include "ps.h"

#include "fit-utils-offdiag.h"
#include "training-set-offdiag.h"

#include "ga.h"
#include "ga_engine.h"
#include "fit_obj.h"
#include "fit_offdiag.h"

#include "dsyev.h"

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

const double ref = -86.15680186351; // two-body interaction energy reference
const double  min_h3o = -76.5300622174*constants::Eh_kcalmol; // hartrees mp2/adz
const double  min_h2o = -76.260909766264*constants::Eh_kcalmol; // hartrees mp2/adz
const double  min_int = -0.137299574*constants::Eh_kcalmol; //hartrees mp2/adz

typedef offdiag::offdiag_v1<degree> model_type;
static model_type model;

hydronium::x1b_hydronium_v1<6> pot_h3o;
h3o_h2o::x2b_h3o_h2o_v1<4> pot_h3o_h2o;

#ifdef RIDGE_REGRESSION
const double alpha = 0.0005;
#endif
const double E_range = 25/constants::Eh_kcalmol; // kcal/mol
    
class FitObj** fit;


static std::vector<offdiag::mol> training_set;
static double*                 ts_weights = 0;

//----------------------------------------------------------------------------//

void diagonalizeH(int nstates, double **H, double &gsenergy) {

  // ** Allocate Evecs and Evals as needed ** //
  double *Evals = new double[nstates];

  double **Evecs = new double*[nstates];
  for (int i=0; i<nstates; ++i) {
    Evecs[i] = new double[nstates];
  }
/*
    for (int i=0; i<nstates; i++) {
        for (int j=0; j<nstates; j++) {
            std::cout<< i << " " << j << " " << H[i][j] << std::endl;
            H[i][j]=1.0;
        }
    }
*/
  // ** Now call the handy diagonalization wrapper ** //
  //printf("Diagonalizing H...\n");
  dsyev(H, nstates, Evals, Evecs);

  // ** Set which eigenvalue and eigenvector is the ground state, store ** //
  // Evals is sorted above such that minimum is element 0, Evecs corresponds

  gsenergy = Evals[0];

  delete [] Evals;
  for (int i=0; i<nstates; ++i) {
    delete [] Evecs[i];
  }
  delete [] Evecs;

}
    
////////////////////////////////////////////////////////////////////////////////

namespace linear {

//----------------------------------------------------------------------------//

static double* A(0);
static double* y(0);
static double* params(0);

//----------------------------------------------------------------------------//

void allocate()
{
    A = new double[training_set.size()*model_type::nparams()];
    y = new double[training_set.size()];

    params = new double[model_type::nparams()];
}

//----------------------------------------------------------------------------//
  
void find_zundel(std::vector<std::string> elements,
                 double *coord,
                 int *fragment,
                 int state,
                 int fragi,
                 int fragj,
                 double *h3o,
                 double *h2o)
{

    std::vector<double> moli(4*3, 0); // hydronium
    std::vector<double> molj(3*3, 0); // water

    int offseti=1;
    int offsetj=1;

    for (int j=0; j<elements.size(); ++j) {
        if (fragment[state*elements.size()+j]==fragi) {
                    
            if (elements[j]=="O") {
                for (int l=0; l<3; ++l) {
                    h3o[0+l] = coord[j*3+l];
                    //moli[0+l] = coord[j*3+l];
                    //                            std::cout << molecule[0+l] << std::endl;}
                    
                }
            }
            else {
                for (int l=0; l<3; ++l){
                    h3o[offseti*3 + l] = coord[j*3+l];
                    //molj[offseti*3 + l] = coord[j*3+l];
                    //                            std::cout << molecule[offset+l] << std::endl;}
                }
                offseti++;
            }
        } else if (fragment[state*elements.size()+j]==fragj) {
                
            if (elements[j]=="O") {
                for (int l=0; l<3; ++l) {
                    h2o[0+l] = coord[j*3+l];
                    //molj[0+l] = coord[j*3+l];
                    //                            std::cout << molecule[0+l] << std::endl;}
                    
                }
            }
            else {
                for (int l=0; l<3; ++l){
                    h2o[offsetj*3 + l] = coord[j*3+l];
                    //molj[offsetj*3 + l] = coord[j*3+l];
                    //                            std::cout << molecule[offset+l] << std::endl;}
                }
                offsetj++;
            }
        }
        
    }

}
    
double compute_coupling(class FitObj** fit,
                        const int n,
                        const int I,
                        const int J)
{
    // Make things shorter hand here
    offdiag::mol ts = training_set[n];
    
    int natoms = ts.natoms;
    double *xx = ts.xyz;
//    char *sym  = mbvb->atom->symbol;
    int *react = ts.reactive;
    double cut_OH  = 2.4;
    double cut_OH2 = cut_OH * cut_OH;
    double inv_cut_OH6 = 1.0 / (cut_OH2 * cut_OH2 * cut_OH2);
    
    double ecouple = 0.0;
    
    // If not shared atom b/w I and J, then coupling is zero
    int shared = 0;
    int iproton = -1; // index of shuttling proton b/w states I and J
    for (int i=0; i<natoms; ++i) {
        if (react[I*natoms + i] && react[J*natoms + i] && ts.elements[i] == "H") {
            // Reactive fragment in I shares at least one H atom with reactive fragment in J
            shared = 1;
            iproton = i;
            break;
        }
    }
    
    if (shared) {
        // Must check that I and J reactive fragments are not on the same oxygen,
        // as it is not possible for self-coupling
        // This can occur with bifurcated waters
        for (int i=0; i<natoms; ++i) {
            if (react[I*natoms + i] && react[J*natoms + i] && ts.elements[i] == "O") {
                // Same oxygen! No coupling allowed.
                shared = 0;
                break;
            }
        }
    }
    
    // ** If passed above, now actually compute the coupling ** //
    if (shared) {
        // Compute interstate coupling repulsion
        double erepItoJ = 0.0;
        for (int i=0; i<natoms; ++i) {
            if (ts.elements[i] == "H" && react[I*natoms + i] && i == iproton) { // reactive H in state I
                for (int j=0; j<natoms; ++j) {
                    if (ts.elements[j] == "O" && react[J*natoms + j]) { // reactive O in state J

                        double dx = xx[3*i]   - xx[3*j];
                        double dy = xx[3*i+1] - xx[3*j+1];
                        double dz = xx[3*i+2] - xx[3*j+2];
                        double dd = dx*dx + dy*dy + dz*dz;
                        if (dd < cut_OH2) {
                            
                            int fragi = ts.fragment[I*natoms+i];
                            int fragj = ts.fragment[I*natoms+j];
                            
                            double h3o[12];
                            double h2o[9];
                            
                            find_zundel(ts.elements,xx,ts.fragment,I,fragi,fragj,h3o,h2o);
                            
                            double xyz[21];
                            std::copy(h3o, h3o + 12, xyz + 0 );
                            std::copy(h2o, h2o + 9 , xyz + 12);
                            
                            ecouple = model.value(xyz)/constants::Eh_kcalmol;
                            //std::cout << "ecouple " << ecouple << std::endl;
                            //std::cout << "zundel" << std::endl;
                            //erepJtoI = pot_offdiag.value(crd)/constants::Eh_kcalmol;
                            

                            //double tmp = Sr(dd, 0.0, cut_OH2);
                            //erepItoJ += tmp*tmp;
                        }
                    }
                }
            }
        }
        // Geometric mean
        //    ecouple = erepItoJ;
        //    ecouple = sqrt ( erepItoJ * erepJtoI );
        //    ecouple = mbvb->math->AA * sqrt(erepItoJ * erepJtoI);
    }
    
    return ecouple;
}

double compute_chisq(const gsl_vector* X, void* unused)
{
    model.set_nonlinear_parameters(X->data);

    if (model.nonlinear_parameters_out_of_range())
        return 1.0e+6;

    PR(model_type::nparams());

    double emin = 1.0e+6;
    for (size_t n = 0; n < training_set.size(); ++n) {
        double mmm[model_type::nparams()];

	model.basis(training_set[n].xyz, mmm);
        y[n] = training_set[n].energy_total_ref
	       - model.basis(training_set[n].xyz, mmm); 
        for (size_t p = 0; p < model_type::nparams(); ++p)
            A[p + n*model_type::nparams()] = mmm[p];
    }

#   ifdef VERBOSE
    std::cout << "=== calling wlsq::solve() ["
              << kit::wlsq::implementation()
              << "] ===" << std::endl;
#   endif

    double chisq;

#   ifdef RIDGE_REGRESSION
    double penaltysq;
    kit::rwlsq::solve(training_set.size(), model_type::nparams(),
                      A, y, ts_weights, alpha, params, chisq, penaltysq);

    std::cout << "<#> chisq = " << chisq
              << " : penaltysq = " << penaltysq
              << std::endl;
#   else
    int rank;
    kit::wlsq::solve(training_set.size(), model_type::nparams(),
                     A, y, ts_weights, params, chisq, rank);
    std::cout << "<#> chisq = " << chisq
              << " : rank = " << rank
              << std::endl;
#   endif

#   ifdef VERBOSE
    std::cout << "\n--> chisq = " << chisq
              << "\n-->  rank = " << rank
              << '\n' << std::endl;
#   endif

    return chisq;
}

//----------------------------------------------------------------------------//

double compute_chisq(class FitObj** fit)
{
    
    int nnlin = model_type::num_nonlinear_params;
    int ncoef = model_type::nparams();
    
    double data[nnlin];
    // Fix Me ...
    std::copy(fit[0]->par, fit[0]->par + nnlin , data);
    // Fix Me ...   
 
    model.set_nonlinear_parameters(data);
    
    if (model.nonlinear_parameters_out_of_range())
        return 1.0e+6;
    
//    PR(model_type::nparams());
    
    double poly[ncoef];
    std::copy(fit[0]->par + nnlin, fit[0]->par + nnlin + ncoef, poly);
    
    //for (size_t l=0; l<ncoef; ++l) std::cout << poly[l] << std::endl;
    model.set(poly);
    double chisq(0);
    
    double emin = 1.0e+6;
    for (size_t n = 0; n < training_set.size(); ++n) {
        //std::cout << std::endl;
        //std::cout << model.value(training_set[n].xyz) << std::endl;
        
        for (int I=0; I<training_set[n].nstates; ++I) {
            for (int J=I+1; J<training_set[n].nstates; ++J) {
                training_set[n].H[I][J] = compute_coupling(fit,n,I,J);
                training_set[n].H[J][I] = training_set[n].H[I][J];
            }
        }
        
        diagonalizeH(training_set[n].nstates,training_set[n].H,training_set[n].energy_total_com);
        
        //model.basis(training_set[n].xyz, poly);
        double tmp = (training_set[n].energy_total_ref - training_set[n].energy_total_com);
//	std::cout << training_set[n].energy_total_ref << " " <<  training_set[n].energy_total_com << std::endl;
//	std::cout << tmp*tmp*ts_weights[n] << std::endl;
        chisq += tmp*tmp*ts_weights[n];

        
    }
        
//    exit(1);
/*
    // X contains the nonlinear parameters
    double alpha_all_atm_types[num_nonlinear_params];
    for(size_t i = 0; i < num_nonlinear_params; ++i){
//        if(X->data[i] < 0)
        if (fit[0]->par[i] < 0)
            return 100000;
//        alpha_all_atm_types[i] = X->data[i];
        alpha_all_atm_types[i] = fit[0]->par[i];
    }

    ttm::smear_ttm4x smr;
    smr.m_aDD_intra_12 = 0.626;
    smr.m_aDD_intra_13 = 0.055;

    double chisq(0);
    for (size_t n = 0; n < training_set.size(); ++n){
        double alpha_this_molec[training_set[n].natm];
        build_atomic_polarizabilities(num_nonlinear_params, alpha_all_atm_types,
                                      training_set[n].atm_type,
                                      alpha_this_molec);

        training_system[n].set_pol(alpha_this_molec);

        ttm::electrostatics elect;
        elect(training_system[n], smr);

        double pol[9];
        elect.polarizability_tensor(training_system[n].nsites, pol);
        chisq += dot_diff_6_vs_9(training_set[n].polarizability, pol);
    }
 

//    printf("[chi]^2 = %lf\n",chisq);
    return chisq;
 */
//    printf("[chi]^2 = %lf\n",chisq);
    return chisq;
}

//----------------------------------------------------------------------------//

double evaluate(int nfit, int igenome, GA_Engine* GA)
{
  double *gene = GA->child[igenome];
  int id=0;

  for(size_t i=0; i<nfit; ++i)
    for(size_t j=0; j<fit[i]->npar; ++j)
      fit[i]->par[j] = gene[id++];

  return compute_chisq(fit);
//  cal_error();
//  return global.error;

}

} // namespace linear

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &(global.nproc));
    MPI_Comm_rank(MPI_COMM_WORLD, &(global.iproc));

    //
    // initialize the model and load training sets
    //

    if (argc < 3) {
        std::cerr << "usage: fit-offdiag-v1 kvs ts1 ..."
                  << std::endl;
        return 0;
    }

    ++argv; --argc;

    if (**argv == 'x') {
	const double x0[] = {
             0.1,
             0.1,
             1.9,
             1.9,
             0.72
        };
//        ++argv;--argc;
//        model.setup(*argv);
        model.set_nonlinear_parameters(x0);
    } else
        model.setup(*argv);

//    pot_h3o.load_netcdf("/Users/oosode/Desktop/offdiag/hydronium.nc");
//    pot_h3o_h2o.load_netcdf("/Users/oosode/Desktop/offdiag/zundel.nc");

#   ifdef RIDGE_REGRESSION
    std::cout << "<> using ridge regression with alpha = "
              << alpha << std::endl;
#   endif

    std::cout << "\n<><><> model type = '" << model_type::name() << "'\n";

    {
        const char fn[] = "fit-offdiag-v1-initial.cdl";
        std::ofstream ofs(fn);
        std::cout << "\n>> dumping initial model as '" << fn << "' >>\n\n";
        model.as_cdl(ofs);
    }

    ++argv;

    try {
        while (--argc != 0) {
            size_t nd = offdiag::load_mols(*argv, training_set);
            std::cout << "'" << *(argv++) << "' : "
                      << nd << " mols" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    //
    // assign weights
    //

    ts_weights = new double[training_set.size()];

    double E_min;
    size_t N_eff;

    offdiag::setup_weights(training_set, E_range, E_min,
                       ts_weights, N_eff, false);

    std::cout << "\n>>   E_min = " << E_min << " kcal/mol"
                 "\n>> E_range = " << E_range << " kcal/mol"
                 "\n\n>> training set size = " << training_set.size()
              << "\n>>    effective size = " << N_eff << '\n'
              << std::endl;
 
    //
    // set things up
    //

    // Read fitting data
    // ** Hard coded values ** //
    int nfit = 1;

    fit = new FitObj* [nfit];
    double *tmp;

    for (size_t i=0; i<nfit; ++i) {
        fit[i] = new FitOffdiag;

        fit[i]->npar = model_type::num_nonlinear_params + model_type::nparams();


        // Allocate memory for parameters
        fit[i]->par = new double [fit[i]->npar];
        fit[i]->min = new double [fit[i]->npar];
        fit[i]->max = new double [fit[i]->npar];
        tmp = new double [fit[i]->npar];     
    }

    if (1) {
        printf("Reading input file %s\n", "tmp.ga.inp");

        FILE *ft = fopen("tmp.ga.inp", "r");

        char line[MAX_LENGTH];

        int l = 0;
        while ( fgets(line, MAX_LENGTH, ft) != NULL ) {
  
            char arg0[MAX_LENGTH];
 
            if ( sscanf(line, "%s", arg0) == 1 ) {
                fit[0]->par[l] = atof(arg0);
  	        tmp[l] = atof(arg0);
                l++;
            }
        }
    }
    // Read min and max for each parameter
    for (size_t i=0; i<nfit; ++i) {
        for(int j=0; j<fit[i]->npar; ++j) {
            if (j<model_type::num_nonlinear_params) {
//                fit[i]->min[j] = -15.0;
//                fit[i]->max[j] =  15.0;
                fit[i]->min[j] =  tmp[j] - 0.25;
                fit[i]->max[j] =  tmp[j] + 0.25;
            } else {
//                fit[i]->min[j] = -300.0;
//                fit[i]->max[j] =  300.0;
                fit[i]->min[j] =  tmp[j] - 0.75;
                fit[i]->max[j] =  tmp[j] + 0.75;
            }
        }
    }
    // **  ** //

    ////////////////////////////////////////////////////////////////////
    if (0) 
    {
    double xi = linear::compute_chisq(fit);
    printf("chisq: %15.10f\n",xi);
    exit(0);
    }
    ////////////////////////////////////////////////////////////////////

    GA_Engine* GA = new GA_Engine();
//    GA->execute();

    printf("Initiate GA fitting ...\n");
    // Set gene
    for (size_t i=0; i<nfit; ++i) GA->ngene += fit[i]->npar;

    GA->allocate_gene();
   
    int id=0;

    for (size_t i=0; i<nfit; ++i) {
      for (size_t j=0; j<fit[i]->npar; ++j) {

        GA->gene_max[id]   = fit[i]->max[j];
        GA->gene_min[id++] = fit[i]->min[j];
      }  
    }    

    GA->begin_init();

    printf("Setup evaluation ...\n");

    for(int i=0; i<GA->nchild; i++) 
    {
      GA->error[i] = linear::evaluate(nfit,i,GA);
    }
    GA->finish_init();

//    timer_reset();
    GA->last_total = 0.0;
    GA->generation = 0;

    printf("Start GA fitting ...\n");

    while(true)
    {
      GA->generation++;

      {
        GA->propagate();
        GA->mutate();
      }
//      timer_click(EVOL);


//      timer_click(COMM);

      for(int i=0; i<GA->nchild; i++) 
      {
        GA->error[i] = linear::evaluate(nfit,i,GA);
//	printf("%10.5f %10.5f\n",GA->error_node[i],nonlinear::evaluate(nfit,i,GA)); 
      }
//      exit(1);
//      timer_click(EVAL);

//    MPI_Reduce(error_node, error, nchild, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//      timer_click(COMM);

//      if(iproc==0) evolute();
      GA->evolute();
//      timer_click(EVOL);

//      if(iproc==0 && generation%every_restart==0)
      if (GA->generation%GA->every_restart==0) GA->output();
//      timer_click(OUTP);

//      if(iproc==0 && (time_stamp - last_total > every_total))
//      if (GA->time_stamp - GA->last_total > GA->every_total)
      {
//        last_total = time_stamp;
        GA->total();
      }

      if(GA->maxgen && GA->generation==GA->maxgen) break;
//      printf("iteration number: %4d\n",GA->generation);

//#if defined(BGQ)
//    if(iproc==global.bg_mem_proc_id) {
//      fprintf(stdout,"\n\nget_memory() at end of GA::execute() generate= %i\n",generation);
//      get_memory();
//    }
//#endif
    }


    std::cout << "Final parameters:" << std::endl;
//    for(size_t i = 0; i < num_nonlinear_params; ++i)
//        std::cout << all_atm_types[i]<< ' ' << s->x->data[i] << std::endl;

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

        double emin = 1.0e+6;
        for (size_t n = 0; n < training_set.size(); ++n) {

            for (int I=0; I<training_set[n].nstates; ++I) {
                for (int J=I+1; J<training_set[n].nstates; ++J) {
                    training_set[n].H[I][J] = linear::compute_coupling(fit,n,I,J);
                    training_set[n].H[J][I] = training_set[n].H[I][J];
                }
            }

            diagonalizeH(training_set[n].nstates,training_set[n].H,training_set[n].energy_total_com);

        }

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


    delete GA;

    for (size_t i=0; i<nfit; ++i) {
	delete [] fit[i]->par;
	delete [] fit[i]->min;
        delete [] fit[i]->max;
//	delete [] fit[i];
    }
    delete [] fit;

    MPI_Finalize();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
