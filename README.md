# BIRRPf
## *WIP*
## Author: Taylor Viti (tviti@hawaii.edu)

Julia bindings to the original BIRRP Fortran subroutines. This really isn't intended to be "useful", moreso as a way of comparing my own implementations to the BIRRP ones.

More info about BIRRP (including the original source) can be obtained at: https://www.whoi.edu/science/AOPE/people/achave/Site/Next1.html

Unfortunately, per the above website, I'm not supposed to be re-distributing the BIRRP source, so the original fortran files need to be placed in the `src/fortran` dir, and then compiled down to a shared like with either (while in `src/fortran`)

`make birrp.dylib` (on OS X), or 

`make birrp.so` (on a *nix).


Users should also be mindful that the `parameters.h` file that came with BIRRP has been configured properly (segfaults are a pretty good sign that it hasn't been).
