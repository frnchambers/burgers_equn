CCCXX  = g++
FORT = gfortran
CXXFLAGS = -Wall -Werror -std=c++11 #-O3
LDFLAGS  = -lgsl -lgslcblas

#################################################################
## -------------------------- Build -------------------------- ##

ROOT_DIR = ../../../src
VPATH = $(ROOT_DIR):$(ROOT_DIR)/FD_derivs:$(ROOT_DIR)/Grid
CXXFLAGS += -I$(ROOT_DIR)

exe = run_grid
exe_bins := $(addsuffix .out, $(exe) )

$(exe_bins): Grid.hpp

## how to build exes
%.out: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -rf *.o $(EXE) data/*.dat
