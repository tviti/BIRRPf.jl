FC = gfortran
FCFLAGS = -std=legacy -ffixed-line-length-0
FCFLAGS += -fbounds-check -march=native
FCFLAGS += -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8
FCFLAGS += -fPIC
LDFLAGS_DYLIB = -dynamiclib
LDFLAGS_SO = -shared

# Useful for debugging (NOTE: BIRRP is likely to throw an FPE, even
# for during "normal" BIRRP behaviour. Go figure...)
# FCFLAGS += -ffpe-trap=invalid,zero,overflow,underflow,denormal

TARGETS = clean  birrp.dylib birrp.so
OBJ = rarfilt.o fft.o math.o utils.o rtpss.o coherence.o zlinpack.o response.o weight.o dataft.o filter.o parameters.o

# These are ALL the object files needed to compile a birrp exec
# OBJ = math.o  utils.o  coherence.o  dataft.o  diagnostic.o fft.o filter.o  rarfilt.o  response.o  rtpss.o    weight.o  zlinpack.o  birrp.o

all:  $(TARGETS)

clean:
	rm -f *.o ./birrp.dylib ./birrp.so

# Compile to a .dylib shared library for OS X usage
# I think this will fail on a *nix system though...
birrp.dylib: $(OBJ)
#	$(FC) $(FCFLAGS) -o parameters.o parameters.h
	$(FC) $(LDFLAGS_DYLIB) -o $@ $(OBJ)

# Compile to a .so shared library for *nix usage
birrp.so: $(OBJ)
	$(FC) $(LDFLAGS_SO) -o $@ $(OBJ)

# 
# Pattern matching rules for Fortran and C files (these replace old style .SUFFIXES rules):
# 
# General Fortran compile:
%.o: %.f 
	$(FC) $(FCFLAGS)    -c -o $@ $^

# %.o: %.f90
# 	$(FC) -std=legacy -ffixed-line-length-72    -c -o $@ $^
