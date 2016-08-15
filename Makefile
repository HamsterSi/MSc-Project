EXE = edmddna

CXX = mpicxx
CC = mpicc
CFLAGS = -pg -O3
LIBS = #-lm

DEP = src/qcprot/qcprot.c
SRC = src/array.cpp src/mpilib.cpp src/tetrad.cpp src/io.cpp src/edmd.cpp src/master.cpp src/worker.cpp src/simulation.cpp src/qcprot/qcprot.c
OBJ1 = $(DEP:.c=.o)
OBJ2 = $(SRC:.cpp=.o) 

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

$(EXE): src/main.cpp $(OBJ2) $(OBJ1)
	$(CXX) $(CFLAGS) $(LIBS) -o $@ $^

clean:
	rm -f src/*.o src/qcprot/qcprot.o $(EXE)
