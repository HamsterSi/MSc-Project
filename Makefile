EXE = main

CXX = mpicxx
CFLAGS = #-g -O3
LIBS = #-lm

SRC = src/simulation.cpp src/tetrad.cpp src/io.cpp src/edmd.cpp src/mpilibrary.cpp src/master.cpp src/worker.cpp 
OBJ = $(SRC:.cpp=.o)

%.o: %.c
	$(CXX) -c $@ $< 
	
$(EXE): src/main.cpp $(OBJ)
	$(CXX) $(CFLAGS) $(LIBS) -o $@ $^ 

clean:
	rm -f src/*.o $(EXE)
