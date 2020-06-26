module AEM_VMD_HMD
include("DigFilters.jl")
using Dierckx, LinearAlgebra

abstract type HField end

mutable struct HFieldDHT <: HField
    thickness   :: Array{Float64, 1}
    pz          :: Array{Complex{Float64}, 1}
    epsc        :: Array{ComplexF64, 1}
    zintfc      :: Array{Float64, 1}
    rTE         :: Array{ComplexF64, 1}
    rTM         :: Array{ComplexF64, 1}
    zRx         :: Float64
    zTx         :: Float64
    rTX         :: Float64
    rRx         :: Float64
    freqs       :: Array{Float64, 1}
    H           :: Array{ComplexF64, 1}
    J0_kernel_h :: Array{ComplexF64, 2}
    J1_kernel_h :: Array{ComplexF64, 2}
    J0_kernel_v :: Array{ComplexF64, 2}
    J1_kernel_v :: Array{ComplexF64, 2}
end

function HFieldDHT(;
      nmax      = 200,
      rTx       = 12.0,
      rRx       = 17.0,
      freqs     = [0.1, 0.3, 0.7],
      zTx       = -35.0,
      zRx       = -37.5
  )
    @assert all(freqs .> 0.)
    thickness = zeros(nmax)
    zintfc    = zeros(nmax)
    pz        = zeros(Complex{Float64}, nmax)
    epsc      = similar(pz)
    rTE       = zeros(length(pz)-1)
    rTM       = similar(rTE)
    H         = zeros(ComplexF64, length(freqs)) # space domain fields
    # define kr domain fields here for lagged convolution
    J0_kernel_h, J1_kernel_h, J0_kernel_v, J1_kernel_v = map(x->zeros(ComplexF64, length(Filter_base), length(freqs)), 1:4)
    HFieldDHT(thickness, pz, epsc, zintfc, rTE, rTM, zRx, zTx, rTx, rRx, freqs, H,
            J0_kernel_h, J1_kernel_h, J0_kernel_v, J1_kernel_v)
end

const mu       = 4*pi*1e-7
const eps0     = 8.854e-12
const iTxLayer = 1
const iRxLayer = 1

function stacks!(F::HField, iTxLayer::Int, nlayers::Int, omega::Float64)

    rTE              = F.rTE
    rTM              = F.rTM
    pz               = view(F.pz, 1:nlayers)
    #The last and first layer thicknesses are infinite
    d                = view(F.thickness, 1:nlayers)
    d[1]             = 1e60
    d[nlayers]       = 1e60

    # Capital R is for a stack
    # Starting from the bottom up, for Rs_down
    Rlowerstack_TE, Rlowerstack_TM = 0. *im, 0. *im
    for k = (nlayers-1):-1:iTxLayer
      Rlowerstack_TE = lowerstack(Rlowerstack_TE, pz, rTE, d, k, omega)
      Rlowerstack_TM = lowerstack(Rlowerstack_TM, pz, rTM, d, k, omega)
    end

return Rlowerstack_TE, Rlowerstack_TM
end

function lowerstack(Rlowerstack::ComplexF64, pz::SubArray{ComplexF64, 1},
                    r::Array{ComplexF64, 1}, d::SubArray{Float64, 1}, k::Int, omega::Float64)
    e_to_the_2iwpznext_dnext = exp(2im*omega*pz[k+1]*d[k+1])
    Rs_d = (r[k] + Rlowerstack * e_to_the_2iwpznext_dnext) /
        (1. + r[k]*Rlowerstack * e_to_the_2iwpznext_dnext)
end

function getCurlyR(Rs_d::ComplexF64, pz::ComplexF64,
                  zR::Float64, z::SubArray{Float64, 1}, iTxLayer::Int, omega::Float64)

    if (zR>=0)
        e_to_the_iwpzzr               = exp( im*omega*pz*zR)
        e_to_the_iwpz_2znext_minus_zr = exp( im*omega*pz*(2*z[iTxLayer+1] - zR))

        finRA = e_to_the_iwpzzr + e_to_the_iwpz_2znext_minus_zr*Rs_d

        finRB = e_to_the_iwpzzr + e_to_the_iwpz_2znext_minus_zr*Rs_d

        finRC = e_to_the_iwpzzr + e_to_the_iwpz_2znext_minus_zr*(-Rs_d)

        finRD = e_to_the_iwpzzr + e_to_the_iwpz_2znext_minus_zr*(-Rs_d)
    else
        e_to_the_iwpz2znext              = exp(im*omega*pz*2*z[iTxLayer+1])
        e_to_the_minus_iwpzzr            = exp(-im*omega*pz*zR)

        finRA = (1. + e_to_the_iwpz2znext*Rs_d) * e_to_the_minus_iwpzzr

        finRB = (-1. + e_to_the_iwpz2znext*Rs_d) * e_to_the_minus_iwpzzr

        finRC = (1. + e_to_the_iwpz2znext*(-Rs_d)) * e_to_the_minus_iwpzzr

        finRD = (-1. + e_to_the_iwpz2znext*(-Rs_d)) * e_to_the_minus_iwpzzr
    end

    return finRA, finRB, finRC, finRD
end

getepsc(rho, omega)      = eps0 + 1im/(rho*omega)
getpz(epsc, krho, omega) = sqrt(mu*epsc - (krho/omega)^2)
ztxorignify(z, zTx)      = z - zTx
makesane(pz::Complex)    = imag(pz)  < 0.0 ? ( pz*=-1.) : pz

function getAEM1DKernelsH!(F::HField, krho::Float64, f::Float64, zz::Array{Float64, 1}, rho::Array{Float64, 1})
    nlayers = length(rho)
    z       = view(F.zintfc, 1:nlayers)
    zRx     = ztxorignify(F.zRx, F.zTx)
    omega   = 2. *pi*f
    epsc    = F.epsc
    pz      = F.pz
    # reflection coefficients (downward) for an intfc: pz vertical slowness
    rTE, rTM = F.rTE, F.rTM

    l = 1
    z[l]       = ztxorignify(zz[l], F.zTx)
    epsc[l]    = getepsc(rho[l], omega)
    pz[l]      = makesane(getpz(epsc[l], krho, omega))
    for intfc in 1:nlayers-1
        l = intfc+1
        z[l]       = ztxorignify(zz[l], F.zTx)
        epsc[l]    = getepsc(rho[l], omega)
        pz[l]      = makesane(getpz(epsc[l], krho, omega))
        rTE[intfc] = (pz[intfc] - pz[intfc+1])/(pz[intfc] + pz[intfc+1])
        rTM[intfc] = (epsc[intfc]*pz[intfc+1] - epsc[intfc+1]*pz[intfc]) /
                     (epsc[intfc]*pz[intfc+1] + epsc[intfc+1]*pz[intfc])
        F.thickness[intfc] = z[intfc+1] - z[intfc]
    end

    # TE and TM modes
    Rs_dTE, Rs_dTM   = stacks!(F, iTxLayer, nlayers, omega)
    curlyRA,         = getCurlyR(Rs_dTE, pz[iTxLayer], zRx, z, iTxLayer, omega)
    gA_TE            = mu/pz[iTxLayer]*curlyRA

    # Kernels according to Loseth, without the bessel functions multiplied
    J0v              = krho^3*gA_TE*1im/(omega*mu*4*pi)

    return J0v
end

function getfield!(F::HFieldDHT, z::Array{Float64, 1}, ρ::Array{Float64, 1})
    for (ifreq, freq) in enumerate(F.freqs)
        for (ikr, kr) in enumerate(Filter_base)
            F.J0_kernel_v[ikr,ifreq] = getAEM1DKernelsH!(F, kr/F.rRx, freq, z, ρ)
        end # kr loop
        F.H[ifreq] = dot(F.J0_kernel_v[:,ifreq], Filter_J0)/F.rRx
    end # freq loop
end

end # module
