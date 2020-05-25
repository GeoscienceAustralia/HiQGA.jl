module CSEM1Dkernels
include("DigFilters.jl")
using Dierckx, LinearAlgebra
export stacks, getCurlyR, getCSEM1DKernelsAnisoHED, getCSEM1DanisoHED

mutable struct RadialEr
  thickness :: Array{Float64, 1}
  pz        :: Array{Complex{Float64}, 1}
  zintfc    :: Array{Float64, 1}
  rTE       :: Array{ComplexF64, 1}
  rTM       :: Array{ComplexF64, 1}
  zRx       :: Array{Float64, 1}
  zTx       :: Array{Float64, 1}
  freqs     :: Array{Float64, 1}
end

const mu      = 4 *pi*1e-7
const eps0    = 8.854e-12

function RadialEr(;
      nmax      = 200,
      zTx       = Nan,
      freqs     = [NaN],
      zRx       = [NaN],
  )
  @assert !isnan(zTx)
  @assert !isnan(zRx[1])
  @assert !isnan(freqs[1])
  thickness = zeros(nmax)
  zintfc  = zeros(nmax)
  pz  = zeros(Complex{Float64}, nmax)
  rTE = zeros(length(pz)-1)
  rTM = similar(rTE)
  RadialEr(thickness, pz, zintfc, rTE, rTM, zRx, zTx, freqs)
end

function stacks(F::RadialEr, iTxLayer::Int, nlayers::Int, omega::Float64)

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
      Rlowerstack_TE = lowerstack(Rlowerstack_TE, pz, rTE, d, k)
      Rlowerstack_TM = lowerstack(Rlowerstack_TM, pz, rTM, d, k)
  end

  Rupperstack_TE, Rupperstack_TM = 0. *im, 0.* im
  # Starting from the top down for Rs_up
  for k   = 2:iTxLayer
      Rupperstack_TE = upperstack(Rupperstack_TE, pz, rTE, d, k)
      Rupperstack_TM = upperstack(Rupperstack_TM, pz, rTM, d, k)
  end

  return Rupperstack_TE, Rlowerstack_TE, Rupperstack_TM, Rlowerstack_TM
end

function lowerstack(Rlowerstack::ComplexF64, pz::Array{ComplexF64, 1}, r::Array{ComplexF64, 1}, d::Array{Float64, 1}, k::Int)
  Rs_d = (r[k] + Rlowerstack*exp(2im*omega*pz[k+1]*d[k+1])) /
        (1. + r[k]*Rlowerstack *
        exp(2im*omega*pz[k+1]*d[k+1]))
end

function upperstack(Rupperstack::ComplexF64, pz::Array{ComplexF64, 1}, r::Array{ComplexF64, 1}, d::Array{Float64, 1}, k::Int)
  Rs_u = (-r[k-1] + Rupperstack*exp(2im*omega*pz[k-1]*d[k-1])) /
        (1. - r[k-1]*Rupperstack *
        exp(2im*omega*pz[k-1]*d[k-1]))
end

function getCurlyR(Rs_u::Array{ComplexF64, 1}, Rs_d::Array{ComplexF64, 1}, pz::Array{ComplexF64, 1},
                  zR::Float64, z::Float64, iTxLayer::Int, omega::Float64)
  d=z[iTxLayer+1]-z[iTxLayer]
  if (zR>=0)
      finRA = (1. + exp(-im*omega*pz*2*z[iTxLayer])*Rs_u) *
              (exp(im*omega*pz*zR) + exp( im*omega*pz*(2*z[iTxLayer+1] - zR))*Rs_d) /
              (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))

      finRB = (1. - exp(-im*omega*pz*2*z[iTxLayer])*Rs_u) *
              (exp(im*omega*pz*zR) + exp( im*omega*pz*(2*z[iTxLayer+1] - zR))*Rs_d) /
              (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))
  else
      finRA = (1. + exp(im*omega*pz*2*z[iTxLayer+1])*Rs_d) *
              (exp(-im*omega*pz*zR) + exp(-im*omega*pz*(2*z[iTxLayer] - zR))*Rs_u) /
              (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))

      finRB = (-1. + exp(1im*omega*pz*2*z[iTxLayer+1])*Rs_d) *
              (exp(-im*omega*pz*zR) + exp(-im*omega*pz*(2*z[iTxLayer] - zR))*Rs_u) /
              (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))
  end

  return finRA, finRB
end

function getCSEM1DKernelsEr(F::RadialEr, krho::Float64, f::Float64, z::Array{Float64, 1}, rho::Array{Float64, 1}, rxno::Int)
  nlayers = length(rho)
  z = z .- F.zTx
  zRx     = F.zRx[rxno]-F.zTx
  omega   = 2. *pi*f
  epsc    = eps0 .+ 1im./(rho*omega)
  pz  = F.pz
  zRxFlag = Array{Int}(undef, length(z))
  zTxFlag = similar(zRxFlag)
  for layer in eachindex(z)
    pz[layer]  = sqrt(mu*epsc[layer] - (krho/omega)^2)
    # wavenubmer sanity check
    imag(pz[layer])  < 0.0 && ( pz[layer]*=-1.)
    #where are the Tx and Rx
    zRxFlag[layer] =  round(Int,z[layer]<zTx)
    zTxFlag[layer] =  round(Int,z[layer]<zRx)
  end
  # reflection coefficients (downward) for an intfc: pz vertical slowness
  rTE = F.rTE
  rTM = F.rTM
  for intfc in 1:nlayers-1
    rTE[intfc] = (pz[intfc] - pz[intfc+1])/(pz[intfc] + pz[intfc+1])
    rTM[intfc] = (epsc[intfc]*pz[intfc+1] - epsc[intfc+1]*pz[intfc]) /
                 (epsc[intfc]*pz[intfc+1] + epsc[intfc+1]*pz[intfc])
    F.thickness[intfc] = z[intfc+1] - z[intfc]
  end
  # Find the layer containing the transmitter:
  iTxLayer = sum(zTxFlag)
  if iTxLayer ==  0
    iTxLayer = length(z) # in the last layer
  end
  # Find the layer containing the receiver:
  iRxLayer = sum(zRxFlag)
  if iRxLayer == 0
    iRxLayer = length(z) # in the last layer
  end
  @assert iTxLayer == iRxLayer "Error, transmitter and receiver need to be located in the same layer! Stopping"
  # TE mode
  Rs_uTE, Rs_dTE, Rs_uTM, Rs_dTM  = stacks(F, iTxLayer, nlayers, omega)
  curlyRA,         = getCurlyR(Rs_uTE, Rs_dTE, pz[iTxLayer], zRx, z, iTxLayer, omega)
  gA_TE            = mu/pz[iTxLayer]*curlyRA

  curlyRA, curlyRB = getCurlyR(Rs_uTM, Rs_dTM, pz[iTxLayer], zRx, z, iTxLayer, omega)
  gA_TM            = pz[iTxLayer]/epsc[iTxLayer]*curlyRA
  gB_TM            = curlyRB

  # Kernels according to Loseth, without the bessel functions multiplied
  # Er from HED and VED
  J0                = -krho*gA_TM/4/pi
  J1                = -(gA_TE - gA_TM)/4/pi
  J1V               = krho^2*gB_TM*1im/(4*pi*omega*epsc[iTxLayer]);
  ErKernels         = [J0; J1; J1V]

  return ErKernels
end

end #module
