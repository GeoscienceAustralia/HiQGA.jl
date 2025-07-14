module CSEM1DEr
include("DigFilters.jl")
using DataInterpolations, LinearAlgebra

abstract type RadialEr end

mutable struct RadialErLagged <: RadialEr
    thickness :: Array{Float64, 1}
    pz        :: Array{Complex{Float64}, 1}
    epsc      :: Array{ComplexF64, 1}
    zintfc    :: Array{Float64, 1}
    rTE       :: Array{ComplexF64, 1}
    rTM       :: Array{ComplexF64, 1}
    zRx       :: Array{Float64, 1}
    zTx       :: Array{Float64, 1}
    rRx       :: Array{Float64, 1}
    freqs     :: Array{Float64, 1}
    Er        :: Array{ComplexF64, 2}
    ErhJ0     :: Array{ComplexF64, 2}
    ErhJ1     :: Array{ComplexF64, 2}
    ErvJ1     :: Array{ComplexF64, 2}
    ErBase    :: Array{ComplexF64, 2}
    lambdaR   :: Array{Float64, 1}
    rR        :: Array{Float64, 1}
    cst       :: Float64
    sit       :: Float64
    csb       :: Float64
    sib       :: Float64
end

const mu      = 4*pi*1e-7
const eps0    = 8.854e-12

function RadialErLagged(;
      nmax      = 200,
      zTx       = [975.],
      rRx       = collect(LinRange(500, 5000, 20)),
      freqs     = [0.1, 0.3, 0.7],
      zRx       = [1000.],
      TxDip     = 0.0,
      RxAzim    = 0.0
  )
    @assert all(freqs .> 0.)
    thickness = zeros(nmax)
    zintfc  = zeros(nmax)
    pz  = zeros(Complex{Float64}, nmax)
    epsc = similar(pz)
    rTE = zeros(length(pz)-1)
    rTM = similar(rTE)
    Er  = zeros(ComplexF64, length(rRx),length(freqs)) # space domain fields
    #lagged convolution, good for all Rx-Tx depths same
    rMin, rMax          = extrema(rRx)
    lamMin              = Filter_base[1]/rMax
    lamMax              = Filter_base[end]/rMin
    # the filters are spaced by a multiplicative factor
    filterSpacing       = Filter_base[2] / Filter_base[1]
    nLambda             = ceil(log(lamMax/lamMin)/log(filterSpacing))+1
    lambdaR             = collect(exp.(log(lamMin):log(filterSpacing):log(lamMin)+log(filterSpacing)*(nLambda-1)))
    nExtra              = nLambda - length(Filter_base)
    # ranges corresponding to lambdaR plus nExtra 1 lags in convolution:
    rR                  = rMax*exp.(-(0:nExtra)*log(filterSpacing))
    # define kr domain fields here for lagged convolution
    ErhJ0, ErhJ1, ErvJ1 = map(x->zeros(ComplexF64, length(lambdaR), length(freqs)), 1:3)
    # needed for spline computations in space domain
    ErBase              = zeros(ComplexF64, length(rR), length(freqs))
    # Dip and azimuth
    # Rx Azimuths and Tx Dips
    cst = cosd(TxDip)  # needed for all HED terms
    sit = sind(TxDip)  # needed for all VED terms
    csb = cosd(RxAzim) # needs to be multiplied into all HED/VED terms inside the loop
    sib = sind(RxAzim) # can be multiplied into all non-VED terms outside the loop
    RadialErLagged(thickness, pz, epsc, zintfc, rTE, rTM, zRx, zTx, rRx, freqs, Er,
            ErhJ0, ErhJ1, ErvJ1, ErBase, lambdaR, rR,
            cst, sit, csb, sib)
end

mutable struct RadialErDHT <: RadialEr
    thickness :: Array{Float64, 1}
    pz        :: Array{Complex{Float64}, 1}
    epsc      :: Array{ComplexF64, 1}
    zintfc    :: Array{Float64, 1}
    rTE       :: Array{ComplexF64, 1}
    rTM       :: Array{ComplexF64, 1}
    zRx       :: Array{Float64, 1}
    zTx       :: Array{Float64, 1}
    rRx       :: Array{Float64, 1}
    freqs     :: Array{Float64, 1}
    Er        :: Array{ComplexF64, 2}
    ErhJ0     :: Array{ComplexF64, 2}
    ErhJ1     :: Array{ComplexF64, 2}
    ErvJ1     :: Array{ComplexF64, 2}
    lambda    :: Array{Float64, 1}
    cst       :: Array{Float64, 1}
    sit       :: Array{Float64, 1}
    csb       :: Array{Float64, 1}
    sib       :: Array{Float64, 1}
end

function RadialErDHT(;
      nmax      = 200,
      rRx       = collect(LinRange(500, 5000, 20)),
      freqs     = [0.1, 0.3, 0.7],
      zTx       = 975*ones(length(rRx)),
      zRx       = 1000*ones(length(rRx)),
      TxDip     = zeros(length(zTx)),
      RxAzim    = zeros(length(rRx))
  )
    @assert all(freqs .> 0.)
    @assert length(zTx) == length(zRx)
    thickness = zeros(nmax)
    zintfc  = zeros(nmax)
    pz  = zeros(Complex{Float64}, nmax)
    epsc = similar(pz)
    rTE = zeros(length(pz)-1)
    rTM = similar(rTE)
    Er  = zeros(ComplexF64, length(rRx),length(freqs)) # space domain fields
    lambda = zeros(length(Filter_base))
    # define kr domain fields here for lagged convolution
    ErhJ0, ErhJ1, ErvJ1 = map(x->zeros(ComplexF64, length(lambda), length(freqs)), 1:3)
    # Dip and azimuth
    # Rx Azimuths and Tx Dips
    cst = cosd.(TxDip)  # needed for all HED terms
    sit = sind.(TxDip)  # needed for all VED terms
    csb = cosd.(RxAzim) # needs to be multiplied into all HED/VED terms inside the loop
    sib = sind.(RxAzim) # can be multiplied into all non-VED terms outside the loop
    RadialErDHT(thickness, pz, epsc, zintfc, rTE, rTM, zRx, zTx, rRx, freqs, Er,
            ErhJ0, ErhJ1, ErvJ1, lambda, cst, sit, csb, sib)
end

function stacks!(F::RadialEr, iTxLayer::Int, nlayers::Int, omega::Float64)

    rTE              = F.rTE
    rTM              = F.rTM
    pz               = view(F.pz, 1:nlayers)
    #The last and first layer thicknesses are infinite
    d                = view(F.thickness, 1:nlayers)
    d[1]             = 1e60
    d[nlayers]       = 1e60

    # Capital R is for a stack
    # Starting from the bottom up, for Rs_down
    Rlowerstack_TE, Rlowerstack_TM = Complex(0., 0.), Complex(0., 0.)
    for k = (nlayers-1):-1:iTxLayer
      Rlowerstack_TE = lowerstack(Rlowerstack_TE, pz, rTE, d, k, omega)
      Rlowerstack_TM = lowerstack(Rlowerstack_TM, pz, rTM, d, k, omega)
    end

    Rupperstack_TE, Rupperstack_TM = Complex(0., 0.), Complex(0., 0.)
    # Starting from the top down for Rs_up
    for k   = 2:iTxLayer
        Rupperstack_TE = upperstack(Rupperstack_TE, pz, rTE, d, k, omega)
        Rupperstack_TM = upperstack(Rupperstack_TM, pz, rTM, d, k, omega)
    end
return Rupperstack_TE, Rlowerstack_TE, Rupperstack_TM, Rlowerstack_TM
end

function lowerstack(Rlowerstack::ComplexF64, pz::SubArray{ComplexF64, 1},
                    r::Array{ComplexF64, 1}, d::SubArray{Float64, 1}, k::Int, omega::Float64)
    a  = Rlowerstack*exp(2im*omega*pz[k+1]*d[k+1])
    Rs_d = (r[k] + a) / (1. + r[k]*a)
end

function upperstack(Rupperstack::ComplexF64, pz::SubArray{ComplexF64, 1},
                    r::Array{ComplexF64, 1}, d::SubArray{Float64, 1}, k::Int, omega::Float64)
    a = Rupperstack*exp(2im*omega*pz[k-1]*d[k-1])
    Rs_u = (-r[k-1] + a) / (1. - r[k-1]*a)
end

function getCurlyR(Rs_u::ComplexF64, Rs_d::ComplexF64, pz::ComplexF64,
                  zR::Float64, z::SubArray{Float64, 1}, iTxLayer::Int, omega::Float64)
    d=z[iTxLayer+1]-z[iTxLayer]
    denom = (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))
    if (zR>=0)
      a = exp(-im*omega*pz*2*z[iTxLayer])
      b = exp(im*omega*pz*zR)
      c = exp(im*omega*pz*(2*z[iTxLayer+1] - zR))
      finRA = (1. + a*Rs_u) * (b + c*Rs_d) / denom
      finRB = (1. - a*Rs_u) * (b + c*Rs_d) / denom
    else
      a = exp(im*omega*pz*2*z[iTxLayer+1])
      b = exp(-im*omega*pz*zR)
      c = exp(-im*omega*pz*(2*z[iTxLayer] - zR))
      finRA = (1. + a*Rs_d) * (b + c*Rs_u) /  denom
      finRB = (-1. + a*Rs_d) * (b + c*Rs_u) / denom
    end

    return finRA, finRB
end

getepsc(rho, omega)      = eps0 + 1im/(rho*omega)
getpz(epsc, krho, omega) = unsafesqrt(mu*epsc - (krho/omega)^2)
checkz(z, zcheck)        = round(Int, z<zcheck)
ztxorignify(z, zTx)      = z - zTx
makesane(pz::Complex)    = imag(pz)  < 0.0 ? ( pz*=-1.) : pz
function unsafesqrt(z::ComplexF64)
    # x = real(z)
    # d = sqrt(x^2+imag(z)^2)
    # ComplexF64(sqrt(0.5(d+x)), sqrt(0.5(d-x)))
    x, y = reim(z)
    ρ = x*x + y*y
    ρ=abs(x)+sqrt(ρ)
    k = -1
    ρ += ρ
    ρ = ldexp(sqrt(ρ),k) #sqrt((abs(z)+abs(x))/2) without over/underflow
    ξ = ρ
    η = y
    η=(η/ρ)/2
    if x<0
        ξ = abs(η)
        η = copysign(ρ,y)
    end
    Complex(ξ,η)
end
function getCSEM1DKernelsEr!(F::RadialEr, krho::Float64, f::Float64, zz::Array{Float64, 1}, rho::Array{Float64, 1},
                            rxno::Int, txno::Int)
    nlayers = length(rho)
    z       = view(F.zintfc, 1:nlayers)
    zRx     = ztxorignify(F.zRx[rxno], F.zTx[txno])
    omega   = 2. *pi*f
    epsc    = F.epsc
    pz      = F.pz
    # reflection coefficients (downward) for an intfc: pz vertical slowness
    rTE, rTM = F.rTE, F.rTM

    l = 1
    z[l]       = ztxorignify(zz[l], F.zTx[txno])
    epsc[l]    = getepsc(rho[l], omega)
    pz[l]      = getpz(epsc[l], krho, omega)
    iTxLayer   =  checkz(z[l], 0)
    iRxLayer   =  checkz(z[l], zRx)
    @inbounds @fastmath for intfc in 1:nlayers-1
        l = intfc+1
        z[l]       = ztxorignify(zz[l], F.zTx[txno])
        epsc[l]    = getepsc(rho[l], omega)
        pz[l]      = getpz(epsc[l], krho, omega)
        iTxLayer  += checkz(z[l], 0)
        iRxLayer  += checkz(z[l], zRx)
        rTE[intfc] = (pz[intfc] - pz[intfc+1])/(pz[intfc] + pz[intfc+1])
        rTM[intfc] = (epsc[intfc]*pz[intfc+1] - epsc[intfc+1]*pz[intfc]) /
                   (epsc[intfc]*pz[intfc+1] + epsc[intfc+1]*pz[intfc])
        F.thickness[intfc] = z[intfc+1] - z[intfc]
    end
    # Find the layer containing the transmitter:
    @assert iTxLayer == iRxLayer "Error, transmitter and receiver need to be located in the same layer! Stopping"
    if iTxLayer ==  0
        iTxLayer = length(z) # in the last layer
        iRxLayer = length(z) # in the last layer
    end
    # TE and TM modes
    Rs_uTE, Rs_dTE, Rs_uTM, Rs_dTM  = stacks!(F, iTxLayer, nlayers, omega)
    curlyRA,         = getCurlyR(Rs_uTE, Rs_dTE, pz[iTxLayer], zRx, z, iTxLayer, omega)
    gA_TE            = mu/pz[iTxLayer]*curlyRA

    curlyRA, curlyRB = getCurlyR(Rs_uTM, Rs_dTM, pz[iTxLayer], zRx, z, iTxLayer, omega)
    gA_TM            = pz[iTxLayer]/epsc[iTxLayer]*curlyRA
    gB_TM            = curlyRB

    # Kernels according to Loseth, without the bessel functions multiplied
    # Er from HED and VED
    ErhJ0              = -krho*gA_TM/4/pi
    ErhJ1              = -(gA_TE - gA_TM)/4/pi
    ErvJ1              = krho^2*gB_TM*1im/(4*pi*omega*epsc[iTxLayer])

    return ErhJ0, ErhJ1, ErvJ1
end

function getfield!(F::RadialErLagged, z::Array{Float64, 1}, ρ::Array{Float64, 1})
    ErhJ0, ErhJ1, ErvJ1, ErBase = F.ErhJ0, F.ErhJ1, F.ErvJ1, F.ErBase
    for (ifreq, freq) in enumerate(F.freqs), (ikr, kr) in enumerate(F.lambdaR)
        ErhJ0[ikr,ifreq],
        ErhJ1[ikr,ifreq],
        ErvJ1[ikr,ifreq] = getCSEM1DKernelsEr!(F, kr, freq, z, ρ, 1, 1)
    end # first loop over freq for evaluations at digital filter specified kr

    for (ifreq, freq) in enumerate(F.freqs)
        for (ir, r) in enumerate(F.rR)
        # Extract results for all rR ranges:
            iShift = (1:length(Filter_base)) .+ ir .-1
            ErBase[ir,ifreq] = F.cst*F.csb*dot(ErhJ0[iShift,ifreq], Filter_J0)/r +
                               F.cst*F.csb*dot(ErhJ1[iShift,ifreq], Filter_J1)/r^2 +
                               F.sit*      dot(ErvJ1[iShift,ifreq], Filter_J1)/r
        end

        # fit spline at computed receiver locations
        # Dierckx needs x to be increasing so end:-1:1
        # Er
        splReal = CubicSpline(real(ErBase[end:-1:1,ifreq]), log10.(F.rR[end:-1:1]),
            extrapolation_left = ExtrapolationType.Linear, extrapolation_right = ExtrapolationType.Linear)
        splImag = CubicSpline(imag(ErBase[end:-1:1,ifreq]), log10.(F.rR[end:-1:1]),
            extrapolation_left = ExtrapolationType.Linear, extrapolation_right = ExtrapolationType.Linear)
        # evaluate spline at required rRx
        F.Er[:,ifreq] = splReal.(log10.(F.rRx)) - 1im*splImag.(log10.(F.rRx))
    end # loop over frequencies
    # done with lagged convolution
end

function getfield!(F::RadialErDHT, z::Array{Float64, 1}, ρ::Array{Float64, 1})
    ErhJ0, ErhJ1, ErvJ1 = F.ErhJ0, F.ErhJ1, F.ErvJ1
    for (ifreq, freq) in enumerate(F.freqs)
        for (ir, r) in enumerate(F.rRx)
            copy!(F.lambda, Filter_base/r)
            for (ikr, kr) in enumerate(F.lambda)
                ErhJ0[ikr,ifreq],
                ErhJ1[ikr,ifreq],
                ErvJ1[ikr,ifreq] = getCSEM1DKernelsEr!(F, kr, freq, z, ρ, 1, ir)
            end # kr loop
            F.Er[ir,ifreq] = F.cst[ir]*F.csb[ir]*dot(ErhJ0[:,ifreq], Filter_J0)/r +
                             F.cst[ir]*F.csb[ir]*dot(ErhJ1[:,ifreq], Filter_J1)/r^2 +
                             F.sit[ir]*          dot(ErvJ1[:,ifreq], Filter_J1)/r
        end # rRx loop
    end # freq loop
end

end # module
