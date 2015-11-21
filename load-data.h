#ifndef LOAD_DATA_H
#define LOAD_DATA_H

#include <cstdlib>
#include <vector>
#include "load-mols.h"
//#include "dsyev.h"

namespace shootpt {

void load_acc(std::vector<mol>& ts, char filename[64]);
void load_basin(std::vector<mol>& ts, char filename[64]);

} // namespace shootpt

#endif // LOAD_DATA_H
