# Precompiler options
CPP_OPTIONS= -DHOST=\"LinuxIFC\"\
             -DMPI -DMPI_BLOCK=32000 \
             -Duse_collective \
             -DCACHE_SIZE=16000 \
             -Davoidalloc \
             -Dvasp6 \
             -Duse_bse_te \
             -Dtbdyn \
             -DVASP2WANNIER90

CPP        = gcc -E -P -C -w $*$(FUFFIX) >$*$(SUFFIX) $(CPP_OPTIONS)

FC         = mpif90 -m64
FCL        = mpif90 -m64

FREE       = -ffree-form -ffree-line-length-none -fallow-argument-mismatch

FFLAGS     = -w

OFLAG      = -O2
OFLAG_IN   = $(OFLAG)
DEBUG      = -O0

LIBDIR     = /usr/lib/x86_64-linux-gnu

LLIBS      = -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LLIBS      += -L/wannier90 -lwannier

#INCS       = -I/usr/include -I/usr/include/mkl -I/usr/include/mkl/fftw
INCS       = -I/usr/include/mkl -I/usr/include/mkl/fftw

OBJECTS    = fftmpiw.o fftmpi_map.o  fftw3d.o  fft3dlib.o

OBJECTS_O1 += fftw3d.o fftmpi.o fftmpiw.o
OBJECTS_O2 += fft3dlib.o

# For what used to be vasp.5.lib
CPP_LIB    = $(CPP)
FC_LIB     = $(FC)
CC_LIB     = gcc
CFLAGS_LIB = -O
FFLAGS_LIB = -O1
FREE_LIB   = $(FREE)

OBJECTS_LIB= linpack_double.o getshmem.o

# For the parser library
CXX_PARS   = g++
#LIBS       += parser
#LLIBS      += -Lparser -lparser -lstdc++
LLIBS      += -lstdc++

# Normally no need to change this
SRCDIR     = ../../src
BINDIR     = ../../bin

