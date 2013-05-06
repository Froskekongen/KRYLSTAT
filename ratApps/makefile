SRC = sncndn.cpp

OBJ = $(SRC:.cpp=.o)

OUT = libRatApps.a

DEP = libarprec.a

# compiler
CXX=g++

# C++ compiler flags
CCFLAGS = -Wall -O3 -fPIC

# include directories
INCLUDES =

.cpp.o:
	$(CXX) $(INCLUDES) $(CCFLAGS) -c $< -o $@

$(OUT): $(OBJ)
	ar x $(DEP)
	ar rcs $(OUT) *.o

clean:
	rm -f *.o $(OUT) Makefile.bak
