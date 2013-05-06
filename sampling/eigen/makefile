SRC=main_ratS.cpp

SRCDIR1=../../lanczos/eigen
SRCDIR2=../../ratApps

DEPS=$(SRCDIR1)/ap.o $(SRCDIR1)/alglibmisc.o $(SRCDIR1)/alglibinternal.o $(SRCDIR1)/linalg.o $(SRCDIR2)/libRatApps.a

OUT=TESTTEST

# compiler
CXX=g++

# C++ compiler flags
CCFLAGS = -O2 -fopenmp

$(OUT):
	$(CXX) $(CCFLAGS) $(SRC) $(SRCDIR1)/ap.o $(SRCDIR1)/alglibmisc.o $(SRCDIR1)/alglibinternal.o $(SRCDIR1)/linalg.o $(SRCDIR2)/libRatApps.a -o $@
	
