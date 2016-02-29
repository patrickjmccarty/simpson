## Tested on CentOS 6.4 ##

```
C_OBJS = *.o

# MAC:
#INCLUDES = -I/usr/include/tcl8.5
#LIBRARIES = -lm -framework Accelerate -framework Tcl
#
# WINDOWS:
#INCLUDES = -IC:/Tcl/include -I../CBLAS/src -IC:/Pthreads/include
#LIBRARIES = -lm *.dll -LC:/Pthreads/lib -lpthreadGC2
#
#LINUX with GNU:
#INCLUDES = -I/usr/include/tcl8.5
#LIBRARIES = -lm -lgslcblas -llapack -L/opt/ActiveTcl-8.5/lib -lActivetcl8.5 -lpthread
#
# LINUX with MKL:
#INCLUDES = -I/opt/intel/mkl/10.2.4.032/include
#LIBRARIES = -L/opt/ActiveTcl-8.5/lib -lActivetcl8.5 -L/opt/intel/mkl/10.2.4.032/lib/em64t -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack -Wl,--end-group -liomp5 -lpthread -lm
#
# LINUX with Atlas and ActiveStateTcl
INCLUDES = -I/opt/ActiveTcl-8.5/include 
LIBRARIES = -lm -L/opt/ActiveTcl-8.5/lib -lActiveTcl8.5 -L/usr/lib64 -lfftw3 -L/usr/lib64/atlas -lcblas -llapack -lpthread



FLAGS = -c -O3 -DNO_NFFT
CC = gcc
# to enable MPI use this
#FLAGS = -c -O3 -DNO_NFFT -DMPI
#CC = mpicc

RM = rm
TAR = tar

simpson: $(C_OBJS)
	$(CC) $(C_OBJS) $(LIBRARIES) -o simpson4 
$(C_OBJS):
	$(CC) $(FLAGS) $(INCLUDES) *.c
clean:
	$(RM) -f *.o 
dist:
	$(TAR) cvzf simpson.tgz *.c *.h simpson.xcodeproj Makefile

```