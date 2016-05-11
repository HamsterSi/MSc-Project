EXE = main

CXX = mpicxx

SRC = io.cpp mpilibrary.cpp parameters.cpp
OBJ = $(SRC:.cpp=.o)

%.o: %.c
	$(CXX) -c $@ $<
	
$(EXE): main.cpp $(OBJ)
	$(CXX) -o $@ $^

clean:
	rm -f *.o
