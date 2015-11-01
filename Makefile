CXX = g++

DEBUG = -g
#DEBUG = -g
CXXFLAGS = -Wall $(DEBUG) -I/opt/local/include -L/opt/local/lib -std=c++11 

#LDFLAGS='-L/opt/local/lib -lnetcdf'

LIBS = -lnetcdf -lgsl -lgslcblas -llapack

#all: fit-x2b-h3o-h2o-v1 eval-pot-h3o-h2o test-x2b-h3o-h2o
all: lmax

lmax: lmax.cpp load-mols.o io-xyz.o  
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o lmax.x lmax.cpp load-mols.o io-xyz.o


load-mols.o: io-xyz.h io-xyz.cpp load-mols.h load-mols.cpp 
	$(CXX) -c $(CXXFLAGS) $(LDFLAGS) io-xyz.cpp load-mols.cpp

clean:
	rm -rf *.o *.exe *.dSYM

.PHONY: all clean
