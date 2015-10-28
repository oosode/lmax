#ifndef LOAD_MOLS_H
#define LOAD_MOLS_H

#include <cstdlib>
#include <vector>
//#include "dsyev.h"

namespace shootpt {

struct mol {

    double *xyz;
    int natoms;
    std::vector<std::string> elements;
};    

size_t load_mols(std::vector<mol>& ts);
double distance(const double* r1, const double* r2);

} // namespace shootpt

#endif // LOAD_MOLS_H
