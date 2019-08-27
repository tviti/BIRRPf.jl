# BIRRPf
## *WIP*
## Author: Taylor Viti (tviti@hawaii.edu)

Julia bindings to the original BIRRP Fortran subroutines.

More info about BIRRP (including the original source) can be obtained at: https://www.whoi.edu/science/AOPE/people/achave/Site/Next1.html

Unfortunately, per the above website, I'm not supposed to be re-distributing the BIRRP source, so the original fortran files need to be placed in the `src/fortran` dir. Users should also be mindful that the `parameters.h` file that came with BIRRP has been configured properly (segfaults are a pretty good sign that it hasn't been).
