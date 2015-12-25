# include <cassert>

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "io-xyz.h"
#include "load-data.h"
#include "constants.h"
#include "load-mols.h"

namespace shootpt {


void load_acc(std::vector<mol>& ts, char filename[64])
{

    std::ifstream ifs(filename);
    std::istringstream iss;
    std::string line;

    if (!ifs) {
        std::ostringstream oss;
        oss << "could not open '" << filename << "' for reading";
        throw std::runtime_error(oss.str());
    }
    
    size_t i = 0;

    size_t lineno(0);
    std::getline(ifs, line);
    while (!ifs.eof()) {

        size_t n;
        std::string result;

        iss.clear();

        iss.str(line);
        iss >> n >> result;
//	std::cout << "mols,acc: " << ts[i].n << " " << n << ts[i].xyz[0] << std::endl;
//	ts[i].n = n;


        if (n==ts[i].n) {
	    std::cout << "hello " << n << " "  << result << std::endl;
	    ts[i].conclusive=1; 
	    if (result.compare("inconclusive") == 0) ts[i].conclusive=0; 
	}

//	std::cout << result << " " << ts[n].n << " " << ts[n].conclusive << std::endl;	
//	std::cout << ts[n].conclusive << " loadacc" << std::endl;
        if (iss.fail()) {
            std::ostringstream oss;
            oss << "unexpected text at line " << lineno << " of the XYZ stream";
            throw std::runtime_error(oss.str());
        }

        std::getline(ifs, line);
        ++lineno;
	++i;
    }

}

void load_basin(std::vector<mol>& ts, char filename[64])
{

//    std::cout << filename << std::endl;

    std::ifstream ifs(filename);
    std::istringstream iss;
    std::string line;

    if (!ifs) {
        std::ostringstream oss;
        oss << "could not open '" << filename << "' for reading";
        throw std::runtime_error(oss.str());
    }

    size_t i = 0;
    size_t lineno(0);
    std::getline(ifs, line);
    while (!ifs.eof()) {

	size_t haf, hab, hbf, hbb;
        std::string result;

        iss.clear();

        iss.str(line);
//	std::cout << line << std::endl;
        iss >> result >> haf >> hab >> hbf >> hbb;

//        std::cout << n << " " <<  result << " " << haf << " " << hab << " " << hbf << " " << hbb << std::endl;

        if (result.compare("Inconclusive") == 0) {
//	    std::cout << "Inconclusive " << ts[n].conclusive << std::endl;
	    if (ts[i].conclusive != 0) break;
	    else {
	        ts[i].haf = haf;
		ts[i].hab = hab;
		ts[i].hbf = hbf;
		ts[i].hbb = hbb;
	    }
        } else if (result.compare("Conclusive") == 0) {
//	    std::cout << "Conclusive " << ts[n].conclusive << std::endl;
	    if (ts[i].conclusive != 1) break;
            else {
//		std::cout<< "good" << std::endl;
                ts[i].haf = haf;
                ts[i].hab = hab;
                ts[i].hbf = hbf;
                ts[i].hbb = hbb;
            }   

	}
	
	ts[i].accepted=0;
	if (haf==1 && hbb==1) ts[i].accepted=1;
	if (hbf==1 && hab==1) ts[i].accepted=1;

        if (iss.fail()) {
//	    std::cout << "HelLL" << std::endl;
            std::ostringstream oss;
            oss << "unexpected text at line " << lineno << " of the XYZ stream";
            throw std::runtime_error(oss.str());
        }

        std::getline(ifs, line);
        ++lineno;
	++i;
    }
    
}

} // namespace shootpt
