SRC=main_eigen_det_app.cpp

#COLPACK Folder
CPDIR=/home/erlenda/libs/ColPack-1.0.9

SRCDIR1=../../lanczos/eigen
SRCDIR2=../../ratApps
SRCDIR3=$(CPDIR)/.libs
# SRCDIR4=$(HOME)/libs/arprec-2.2.7/src

OUT=TESTTEST

# compiler
CXX=g++


# C++ compiler flags
CCFLAGS = -O2 -fopenmp -march=native -mtune=native

#Header includes. Needs eigen + COLPACK
INCLUDE=-I$(HOME)/libs/eigen -I$(CPDIR)/GraphColoring -I$(CPDIR)/Utilities -I$(CPDIR)/Main -I$(CPDIR)/BipartiteGraphBicoloring -I$(CPDIR)/BipartiteGraphPartialColoring -I$(CPDIR)/Recovery


#$(OUT):
#	$(CXX) $(CCFLAGS) $(SRC) $(SRCDIR1)/ap.o $(SRCDIR1)/alglibmisc.o $(SRCDIR1)/alglibinternal.o $(SRCDIR1)/#linalg.o $(SRCDIR2)/libRatApps.a $(SRCDIR3)/libColPack.a -o $@

#$(OUT):
#	$(CXX) $(CCFLAGS) $(INCLUDE) $(SRC) $(SRCDIR1)/ap.o $(SRCDIR1)/alglibmisc.o $(SRCDIR1)/alglibinternal.o $(SRCDIR1)/linalg.o $(SRCDIR3)/libColPack.a -o $@

$(OUT):
	$(CXX) $(CCFLAGS) $(INCLUDE) $(SRC) $(SRCDIR1)/ap.o $(SRCDIR1)/alglibmisc.o $(SRCDIR1)/alglibinternal.o $(SRCDIR1)/linalg.o $(SRCDIR3)/libColPack.a -o $@
