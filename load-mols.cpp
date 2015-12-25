# include <cassert>

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "io-xyz.h"
#include "load-mols.h"
#include "constants.h"

namespace shootpt {


size_t load_mols(std::vector<mol>& ts)
{

    // while loop over all relevant shootpt files (PDB)
//    while ()
    size_t nshoots = 0;

    for (int i=1; i<=46; i++)
    {

//      std::cout << i << std::endl;
      char filename[64];
      sprintf(filename, "shootpt/%d.shootpt.pdb",i);
//      const char* filename = "115.shootpt.pdb";
//      assert(filename);
    
      std::ifstream ifs(filename);
      if (!ifs) {
          std::ostringstream oss;
          oss << "could not open '" << filename << "' for reading";
          throw std::runtime_error(oss.str());
      }
    
      size_t nmols(0);
    
      std::string comment;
      std::vector<std::string> elements;
      std::vector<double> coord;
      int natoms;
    
      std::istringstream iss;
    
      kit::io::load_pdb(ifs, natoms, comment, elements, coord);
      if (elements.empty())
          break;

      iss.clear();
      iss.str(comment);
        
//          double energy_total_ref;
//          double energy_twobody;
//          double energy_onebody[2];
          
//          iss >> energy_total_ref;
        
      if (iss.fail()) {
          std::ostringstream oss;
          oss << "'" << filename << "' : configuration #"
              << (nmols + 1) << " : unexpected text '";
              //<< iss << "' instead of total/twobody/onebody[2] energies";
          throw std::runtime_error(oss.str());
      }
       


      mol m;
      m.n = i;
      m.natoms = natoms;
      m.xyz = new double[natoms*3];
//          m.energy_total_ref = energy_total_ref;

      std::copy(coord.begin(), coord.end(), m.xyz);
      m.elements = elements;

      ts.push_back(m);
//      std::cout << "hello" << std::endl;
//      std::cout << m.n << std::endl;
        
      ++nshoots;
    }

    return nshoots;
}

} // namespace shootpt
