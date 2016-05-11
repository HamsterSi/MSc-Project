EXE = main

CXX = mpicxx
CFLAGS = -g -O3
LIBS = -lm

DEP = fortran_matfit/matfit.o fortran_matfit/sum.o pool/pool.o
SRC = tetrad.cpp io.cpp edmd.cpp mpilibrary.cpp simulation.cpp
OBJ = $(SRC:.cpp=.o)

%.o: %.c
	#$(CXX) -c $@ $< 
	$(CXX) -c tetrad.o tetrad.cpp fortran_matfit/matfit.o fortran_matfit/sum.o
	$(CXX) -c io.o io.cpp 
	$(CXX) -c parameters.o parameters.cpp 
	$(CXX) -c edmd.o edmd.cpp
	$(CXX) -c master.o master.cpp pool/pool.o
	$(CXX) -c worker.o worker.cpp pool/pool.o
	$(CXX) -c simulation.o simulation.cpp pool/pool.o
	
$(EXE): main.cpp $(OBJ)
	$(CXX) $(CFLAGS) $(LIBS) -o $@ $^ $(DEP)

clean:
	rm -f *.o $(EXE)
