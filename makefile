##############################################################
#
#		Makefile
#		
#		To compile:
#   random number generator (.o file): make ran
#		program to compute traces:  make main
#		program to compute observables: make anal
#   clean: make clean		
#
###############################################################
AT := LOCAL

ifeq ($(AT),LOCAL)
LIBDIR := ~/lattice/libs/
INCLUDE_MKL := ./include
LIBS := -llapack -lblas
MKLFLAGS := -L$(LIBDIR)
INCLUDE := ./include
endif

ifeq ($(AT),EMMY)
LIBDIR := /installadmin/software/intel/composer_xe_2013.3.163/mkl/lib/intel64/
INCLUDE_MKL := /installadmin/software/intel/composer_xe_2013.3.163/mkl/include/
LIBS := $(LIBDIR)libmkl_lapack95_lp64.a $(LIBDIR)libmkl_blas95_lp64.a
MKLFLAGS := -L$(LIBDIR) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5 -lm
INCLUDE := ./include
endif

ifeq ($(AT),FRODO)
LIBDIR := /installadmin/software/intel/composer_xe_2013.3.163/mkl/lib/intel64/
INCLUDE_MKL := /installadmin/software/intel/composer_xe_2013.3.163/mkl/include/
LIBS := $(LIBDIR)libmkl_lapack95_lp64.a $(LIBDIR)libmkl_blas95_lp64.a
MKLFLAGS := -L$(LIBDIR) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5 -lm
INCLUDE := ./include
endif

CFLAGS := -O3 -Wall -I$(INCLUDE) -I$(INCLUDE_MKL)

ifneq ($(AT),LOCAL)
CFLAGS += -DMKL
endif

ran : ranlxd.c
	@echo "compiling for" $(AT)
	gcc $(CFLAGS) -c ranlxd.c

main: ranlxd.o main.cpp
	@echo "compiling for" $(AT)	
	mpicxx $(CFLAGS) -o ./bin/tr.x main.cpp ranlxd.o $(MKLFLAGS) $(LIBS)

debug: ranlxd.o main.cpp
	@echo "compiling for" $(AT)
	mpicxx $(CFLAGS) -DDEBUG -o ./bin/tr.x main.cpp ranlxd.o $(MKLFLAGS) $(LIBS)

anal: analysis.cpp
	@echo "compiling for" $(AT)
	g++ $(CFLAGS) -o ./bin/anal.x analysis.cpp $(MKLFLAGS) $(LIBS)

free: free_case.cpp
	@echo "compiling for" $(AT)
	g++ -std=c++0x -O3 -o ./bin/free.x free_case.cpp

clean:
	@echo "compiling for" $(AT)
	rm -v *.o ./bin/*.x
