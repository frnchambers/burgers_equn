CXX = g++
CXXFLAGS = -Wall -Werror -std=c++11 #-O3
LDFLAGS = 

EXT_DIR=./src
VPATH = $(EXT_DIR)/Interpolation:$(EXT_DIR)/Grid
CXXFLAGS += -I$(EXT_DIR)

## how to build exes
%.out: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJECTS) $(LDFLAGS)

## how to build objects
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<



clean:
	rm -rf *.o *.out data/*.dat
