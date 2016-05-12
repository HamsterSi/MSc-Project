EXE = edmd

CXX = mpicxx
CFLAGS = -g -O3
LIBS = -lm

DEP = src/fortran_matfit/matfit.o src/fortran_matfit/sum.o src/pool/pool.o
SRC = src/tetrad.cpp src/io.cpp src/edmd.cpp src/mpilibrary.cpp src/master.cpp src/worker.cpp src/simulation.cpp
OBJ = $(SRC:.cpp=.o)

%.o: %.c
	$(CXX) -c $@ $< 
	#$(CXX) -c tetrad.o tetrad.cpp fortran_matfit/matfit.o fortran_matfit/sum.o
	#$(CXX) -c io.o io.cpp 
	#$(CXX) -c parameters.o parameters.cpp 
	#$(CXX) -c edmd.o edmd.cpp
	#$(CXX) -c master.o master.cpp pool/pool.o
	#$(CXX) -c worker.o worker.cpp pool/pool.o
	#$(CXX) -c simulation.o simulation.cpp pool/pool.o
	
$(EXE): src/main.cpp $(OBJ)
	$(CXX) $(CFLAGS) $(LIBS) -o $@ $^ 

clean:
	rm -f src/*.o $(EXE)
