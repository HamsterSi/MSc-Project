EXE = edmddna

CXX = mpicxx
CFLAGS = #-g -O3
LIBS = #-lm

DEP = src/qcprot/qcprot.o
SRC = src/simulation.cpp src/tetrad.cpp src/io.cpp src/edmd.cpp src/mpilibrary.cpp src/master.cpp src/worker.cpp 
OBJ = $(SRC:.cpp=.o)

%.o: %.c 
	#$(CXX) -c $@ $< 
	$(CXX) -c src/simulation.o src/simulation.cpp
	$(CXX) -c src/tetrad.o src/tetrad.cpp
	$(CXX) -c src/io.o src/io.cpp
	$(CXX) -c src/edmd.o src/edmd.cpp 
	$(CXX) -c src/mpilibrary.o src/mpilibrary.cpp
	$(CXX) -c src/master.o src/master.cpp
	$(CXX) -c src/worker.o src/worker.cpp

$(EXE): src/main.cpp $(OBJ) 
	$(CXX) $(CFLAGS) $(LIBS) -o $@ $^

clean:
	rm -f src/*.o $(EXE)
