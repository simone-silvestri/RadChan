

SRC  = ./src/
SDNS = ./src/dns/
SMC  = ./src/montecarlo/
SPAR = ./src/param/
OBJ = ./obj/
MOD = ./mod/

DECOMP_OPTIONS= -DDOUBLE_PREC


FLAGS = -c
ifeq ($(DBGc),-g)
F77 = ftn -132 -r8 -mcmodel=large
else 
F77 = mpif90 -ffixed-line-length-none -fdefault-real-8 -O2 -mcmodel=large  -J $(MOD)
endif

NVCC = nvcc -Xcompiler -mno-float128
FLAG1 = -arch 'compute_70' -code 'sm_70' #-maxrregcount 32
PROF = --ptxas-options=-v  -lineinfo
FLAG2 = --use_fast_math
MAT = -ftz=true -prec-div=false
EXEC = g++ -mcmodel=large 

# linking libraries and includes for heterogeneous cuda-fortran programming
#CUDA = /opt/cuda-8.0
CUDA = /cineca/prod/opt/compilers/cuda/10.1/none/
INC = -I$(CUDA)/include
LIB = -L$(CUDA)/lib64 -lc -lstdc++ -lcuda -lcurand -lcudart

DECOMP = $(HOME)/2decomp_fft
INCD   = -I$(DECOMP)/include
LIBD   = -L$(DECOMP)/lib -l2decomp_fft


PROGRAM = RadChan

OBJS = $(OBJ)channel.o $(OBJ)vfft.o $(OBJ)fileio.o $(OBJ)math.o $(OBJ)ccf.o \
$(OBJ)mk_grid.o $(OBJ)fixvar.o $(OBJ)spline.o $(OBJ)mc_gpu_noint.o $(OBJ)read_tables.o \
$(OBJ)memory_copy_noint.o $(OBJ)routine.o


all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F77) $(OBJS) $(INC) $(LIB) -o $(PROGRAM) $(INCD) $(LIBD)

$(OBJ)spline.o: $(SMC)spline.cpp
	$(EXEC) $(DBGc) $(FLAGS) $(OPT) $(SMC)spline.cpp -o $(OBJ)spline.o
$(OBJ)mc_gpu_noint.o: $(SMC)mc_gpu_noint.cu $(SPAR)NarrowBand.h $(SPAR)param.h
	$(NVCC) $(DBG) $(OPT) $(MAT) $(PROF) $(PTX) $(FLAG2) $(CACHE)  $(FLAGS) $(OPT) $(SMC)mc_gpu_noint.cu $(FLAG1) -o $(OBJ)mc_gpu_noint.o
$(OBJ)memory_copy_noint.o: $(SMC)memory_copy_noint.cu
	$(NVCC) $(DBG) $(OPT) $(MAT) $(PROF) $(PTX) $(FLAG2) $(CACHE)  $(FLAGS) $(OPT) $(SMC)memory_copy_noint.cu $(FLAG1) -o $(OBJ)memory_copy_noint.o
$(OBJ)read_tables.o: $(SMC)read_tables.cpp
	$(EXEC) $(DBGc) $(FLAGS) $(OPT) $(SMC)read_tables.cpp -o $(OBJ)read_tables.o
$(OBJ)fixvar.o : $(SDNS)fixvar.f $(SPAR)param.txt $(SDNS)common.txt
	$(F77) $(DBGc) $(FLAGS)  $(SDNS)fixvar.f $(INCD) $(LIBD) -o $(OBJ)fixvar.o
$(OBJ)routine.o : $(SDNS)routine.f90 $(SPAR)param.txt 
	$(F77) $(DBGc) $(FLAGS)  $(SDNS)routine.f90 $(INCD) $(LIBD) -o $(OBJ)routine.o
$(OBJ)channel.o : $(SDNS)channel.f $(SDNS)param.txt $(SDNS)common.txt
	$(F77) $(DBGc) $(FLAGS) $(SDNS)channel.f $(INCD) $(LIBD) -o $(OBJ)channel.o
$(OBJ)fileio.o : $(SDNS)fileio.f $(SPAR)param.txt $(SDNS)common.txt
	$(F77) $(DBGc) $(FLAGS)  $(SDNS)fileio.f $(INCD) $(LIBD) -o $(OBJ)fileio.o
$(OBJ)math.o : $(SDNS)math.f $(SPAR)param.txt $(SDNS)common.txt 
	$(F77) $(DBGc) $(FLAGS)  $(SDNS)math.f $(INCD) $(LIBD) -o $(OBJ)math.o
$(OBJ)ccf.o : $(SDNS)ccf.f $(SPAR)param.txt $(SDNS)common.txt 
	$(F77) $(DBGc) $(FLAGS)  $(SDNS)ccf.f $(INCD) $(LIBD) -o $(OBJ)ccf.o
$(OBJ)mk_grid.o : $(SDNS)mk_grid.f $(SPAR)param.txt $(SDNS)common.txt
	$(F77) $(DBGc) $(FLAGS)  $(SDNS)mk_grid.f $(INCD) $(LIBD) -o $(OBJ)mk_grid.o
$(OBJ)vfft.o : $(SDNS)vfft.f 
	$(F77) $(DBGc) $(FLAGS) $(SDNS)vfft.f -o $(OBJ)vfft.o


clean:
	$(RM) a.out core $(MOD)*.mod $(OBJ)*.o $(PROGRAM)

