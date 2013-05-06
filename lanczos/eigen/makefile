OBJ_FILES = alglib/linalg.o alglib/ap.o alglib/alglibinternal.o alglib/alglibmisc.o main_lancz.o
LD_FLAGS = ...
CC_FLAGS = -O2 -fopenmp

DEST = main_lancz

$DEST: $(OBJ_FILES)
	g++ $< -o $@

%.o: %cpp
	g++ $CC_FLAGS -o -c $@ $<
