# include <cassert>

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "x2b-h3o-h2o-v1.h"
#include "x1b-hydronium-v1.h"
#include "ps.h"
#include "x2b-v9x.h"
//#include "dsyev.h"

#include "io-xyz.h"
#include "training-set-offdiag.h"
#include "constants.h"
//#include "xyz-water-utils.h"

namespace offdiag {

int MAX_STATES = 10;
int MAX_HOPS   = 3;

const double ref = -86.15680186351; // two-body interaction energy reference
const double  min_h3o_Eh = -76.5300622174;
const double  min_h2o_Eh = -76.260909766264;
const double  min_h3o = -76.5300622174*constants::Eh_kcalmol; // hartrees mp2/adz
const double  min_h2o = -76.260909766264*constants::Eh_kcalmol; // hartrees mp2/adz
const double  min_int = -0.137299574*constants::Eh_kcalmol; //hartrees mp2/adz

hydronium::x1b_hydronium_v1<6> pot_h3o;
h3o_h2o::x2b_h3o_h2o_v1<4> pot_h3o_h2o;


size_t load_mols(std::vector<mol>& ts)
{

    // while loop over all relevant shootpt files (PDB)
    while ()
    {
      assert(filename);
    
      std::ifstream ifs(filename);
      if (!ifs) {
          std::ostringstream oss;
          oss << "could not open '" << filename << "' for reading";
          throw std::runtime_error(oss.str());
      }
    
      hydronium::x1b_hydronium_v1<6> pot_h3o;
  //    h3o_h2o::x2b_h3o_h2o_v1<4> pot_h3o_h2o;
    
      pot_h3o.load_netcdf("/Users/oosode/Desktop/offdiag/hydronium.nc");
  //    pot_h3o_h2o.load_netcdf("zundel.nc");
    
      size_t nmols(0);
    
      std::string comment;
      std::vector<std::string> elements;
      std::vector<double> coord;
      int natoms;
    
      std::istringstream iss;
    
      while (!ifs.eof()) {
          kit::io::load_xyz(ifs, natoms, comment, elements, coord);
          if (elements.empty())
              break;

          iss.clear();
          iss.str(comment);
        
          double energy_total_ref;
          double energy_twobody;
          double energy_onebody[2];
          
          iss >> energy_total_ref;
        
          if (iss.fail()) {
              std::ostringstream oss;
              oss << "'" << filename << "' : configuration #"
              << (nmols + 1) << " : unexpected text '";
              //<< iss << "' instead of total/twobody/onebody[2] energies";
              throw std::runtime_error(oss.str());
          }
       


          mol m;
  	  m.natoms = natoms;
          m.xyz = new double[natoms*3];
          m.energy_total_ref = energy_total_ref;

          std::copy(coord.begin(), coord.end(), m.xyz);
          m.elements = elements;

          m.fragment = new int[MAX_STATES*natoms];
          m.reactive = new int[MAX_STATES*natoms];
        
          find_fragments(elements,coord,m.fragment,m.nfragments,m.ireactive);
          state_search(elements,coord,m.fragment,m.reactive,m.nfragments,m.ireactive,m.nstates);


          // Allocate H based on number of states for this step
          m.H = new double*[m.nstates];
          for (int i=0; i<m.nstates; ++i) {
            m.H[i] = new double[m.nstates];
          }

          populateH(elements,coord,m.fragment,m.nfragments,m.nstates,m.H);

          ts.push_back(m);

        
          ++nshoots;
      }

    }  
    
    return nshoots;
}

} // namespace offdiag
