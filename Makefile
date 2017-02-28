
CXX = g++
#CXXFLAGS = -Wall -Werror -ansi -DNDEBUG -mavx -std=c++14 #-O3
CXXFLAGS =  -ansi -DNDEBUG -mavx -Wall -Werror -std=c++14 $(FLAGS)
LDFLAGS = -lboost_program_options -llapack

#CXXFLAGS += -I/home/frank/bin/blaze-3.1/:/home/frank/Documents/burgers_equn/src/
CXXFLAGS += -I./src/ -I/home/frank/bin/blaze-3.1

## how to build exes
%.out: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -rf *.o *.out data/*.dat
