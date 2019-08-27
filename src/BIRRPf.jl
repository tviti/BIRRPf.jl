#=
Interfaces to Fortran subroutines taken from the original BIRRP
source.
=#

module BIRRPf

import Libdl

# These should correspond to the types that BIRRP was compiled with.
fortfloat = Float64
fortcomplex = Complex{Float64}
fortint = Int64

# "locate" and load the shared library
const _BIRRP_fn = Libdl.find_library(["birrp"], [(@__DIR__)*"/fortran"])
const _birrp_dl = Libdl.dlopen(_BIRRP_fn)

# Grab pointers to all of the functions that we are going to expose
const _rarfilt_sym = Libdl.dlsym(_birrp_dl, :rarfilt_)
const _dpqr79_sym = Libdl.dlsym(_birrp_dl, :dpqr79_)
const _sft_sym = Libdl.dlsym(_birrp_dl, :sft_)
const _prewhiten_sym = Libdl.dlsym(_birrp_dl, :prewhiten_)


"""
Allocate some space for computations that need it (e.g. rarfilt).
"""
function _init_workspace(npts::Int64, nfft::Int64)
    work = zeros(fortfloat, 18*nfft+14*npts+7*ceil(Int64, npts/nfft)+32)

    return work
end


"""
Robustly compute ar filter coefficients for prewhitening a vector of data.
"""
function rarfilt(rdata::Vector{fortfloat}, nar::fortint, nfft::fortint)
    ar = zeros(fortfloat, (nar, ))

    work = _init_workspace(size(rdata, 1), nfft)
    
    rarfilt!(rdata, nar, nfft, ar, work)

    return ar
end


"""
Robustly compute ar filter coefficients for multiple series of data.
"""
function rarfilt(rdata::Array{fortfloat, 2}, nar::fortint, nfft::fortint)
    N_ch = size(rdata, 2)

    ar = zeros(fortfloat, (nar, N_ch))
    ar_tmp = zeros(fortfloat, (nar, ))

    work = _init_workspace(size(rdata, 1), nfft)

    for i in 1:N_ch
        ar[:, i] = rarfilt!(rdata[:, i], nar, nfft, ar_tmp, work)
    end

    return ar
end


function rarfilt!(rdata::AbstractVector{fortfloat},
                  nar::fortint,
                  nfft::fortint,
                  ar::AbstractVector{fortfloat},
                  work::Vector{fortfloat})
     npts = size(rdata, 1)
     istat = Ref{fortint}(0)

     ccall(_rarfilt_sym,
          Cvoid,
          (Ref{fortint},  # npts
           Ref{fortint},  # nfft
           Ref{fortint},  # nar
           Ref{fortfloat},  # rdata
           Ref{fortfloat},  # work,
           Ref{fortfloat},  # ar,
           Ref{fortint}),  # istat
          npts,
          nfft,
          nar,
          rdata,
          work,
          ar,
          istat)

     if istat.x != 0
        error("Error in call to rarfilt! istat = ", istat)
    end

    return ar
end


"""
Compute the zeros of a polynomial defined by its coefficients
"""
function dpqr79(coeff::Vector{fortfloat})
    ndeg = size(coeff, 1) - 1
    root = zeros(fortcomplex, (ndeg, ))
    work = zeros(fortfloat, (ndeg*(ndeg+2) + ndeg, ))

    dpqr79!(coeff, root, work)

    return root
end


function dpqr79!(coeff::AbstractVector{fortfloat},
                 root::AbstractVector{fortcomplex},
                 work::AbstractVector{fortfloat})

    ndeg = size(coeff, 1) - 1
    ierr = Ref{fortint}(0)

    ccall(_dpqr79_sym,
          Cvoid,
          (Ref{fortint},  # ndeg
           Ref{fortfloat},  # coeff
           Ref{fortcomplex},  # root
           Ref{fortint},  # ierr
           Ref{fortfloat}),  # work
          ndeg,
          coeff,
          root,
          ierr,
          work)

    if ierr.x != 0
        error("ierr = ", ierr.x)
    end

    return root
end


function sft(x::AbstractVector{fortfloat}, om::fortfloat)
    ct = Ref{fortfloat}(0.0)
    st = Ref{fortfloat}(0.0)

    sft!(x, om, ct, st)

    return ct.x + im*st.x
end


function sft!(x::AbstractVector{fortfloat},
              om::fortfloat,
              ct::Ref{fortfloat},
              st::Ref{fortfloat})

    n = length(x)

    ccall(_sft_sym,
          Cvoid,
          (Ref{fortfloat},  # x
           Ref{fortint},  # n
           Ref{fortfloat},  # om
           Ref{fortfloat},  # ct
           Ref{fortfloat}),  # st
          x,
          n,
          om,
          ct,
          st)

    return
end


function prewhiten(rdata::AbstractVector{fortfloat},
                   ar::AbstractVector{fortfloat})

    pdata = zeros(fortfloat, (length(rdata), ))

    prewhiten!(rdata, ar, pdata)

    return pdata
end


function prewhiten!(rdata::AbstractVector{fortfloat},
                    ar::AbstractVector{fortfloat},
                    pdata::AbstractVector{fortfloat})

    if length(pdata) != length(rdata)
        error("unequal size input and output arrays")
    end

    npts = length(rdata)
    nar = length(ar)

    ccall(_prewhiten_sym,
          Cvoid,
          (Ref{fortint},  # npts
           Ref{fortint},  # nar
           Ref{fortfloat},  # rdata
           Ref{fortfloat},  # ar
           Ref{fortfloat}),  # pdata
          npts,
          nar,
          rdata,
          ar,
          pdata)

    return pdata
end


end  # module BIRRPf
