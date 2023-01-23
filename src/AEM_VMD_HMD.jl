module AEM_VMD_HMD
include("DigFilters.jl")
include("TDcossin.jl")
using LinearAlgebra, SpecialFunctions, FastGaussQuadrature, DataInterpolations
export getfreqsperdec, gettimesperdec
abstract type HField end

mutable struct HFieldDHT <: HField
    thickness       :: Array{Float64, 1}
    pz              :: Array{Complex{Float64}, 1}
    ϵᵢ              :: Array{ComplexF64, 1}
    zintfc          :: Array{Float64, 1}
    rTE             :: Array{ComplexF64, 1}
    rTM             :: Array{ComplexF64, 1}
    zRx             :: Float64
    zTx             :: Float64
    rTx             :: Union{Float64, Nothing}
    rRx             :: Float64
    freqs           :: Array{Float64, 1}
    times           :: Array{Float64, 1}
    ramp            :: Array{Float64, 2}
    log10ω          :: Array{Float64, 1}
    interptimes     :: Array{Float64, 1}
    HFD_z           :: Array{ComplexF64, 1}
    HFD_r           :: Array{ComplexF64, 1}
    HFD_az          :: Array{ComplexF64, 1}
    HFD_z_interp    :: Array{ComplexF64, 1}
    HFD_r_interp    :: Array{ComplexF64, 1}
    HFD_az_interp   :: Array{ComplexF64, 1}
    HTD_z_interp    :: Array{Float64, 1}
    HTD_r_interp    :: Array{Float64, 1}
    HTD_az_interp   :: Array{ComplexF64, 1}
    dBzdt           :: Array{Float64, 1}
    dBrdt           :: Array{Float64, 1}
    dBazdt          :: Array{Float64, 1}
    J0_kernel_h     :: Array{ComplexF64, 2}
    J1_kernel_h     :: Array{ComplexF64, 2}
    J0_kernel_v     :: Array{ComplexF64, 2}
    J1_kernel_v     :: Array{ComplexF64, 2}
    J01kernelhold   :: Vector
    lowpassfcs      :: Array{Float64, 1}
    quadnodes       :: Array{Float64, 1}
    quadweights     :: Array{Float64, 1}
    interplog10ω    :: Array{Float64, 2}
    ω               :: Array{Float64, 2}
    Hsc             :: Array{ComplexF64, 2}
    rxwithinloop    :: Bool
    provideddt      :: Bool
    doconvramp      :: Bool
    useprimary      :: Float64
    nkᵣeval         :: Int
    interpkᵣ        :: Array{Float64,1}
    log10interpkᵣ   :: Array{Float64,1}
    log10Filter_base   :: Array{Float64,1}
    getradialH      :: Bool
    getazimH        :: Bool
    calcjacobian    :: Bool
    Jtemp           :: Vector
    b               :: Vector # Another Jtemp
    derivmatrix     :: Array{ComplexF64, 2}
    HFD_z_J         :: Array{ComplexF64, 2}
    HTD_z_J_interp  :: Array{Float64, 2}
    dBzdt_J         :: Array{Float64, 2}
    HFD_az_J
    HTD_az_J_interp     
    dBazdt_J
    HFD_r_J
    HTD_r_J_interp     
    dBrdt_J    
end

function HFieldDHT(;
      nmax      = 80,
      rTx       = nothing,
      rRx       = 17.0,
      freqs     = [],
      zTx       = -35.0,
      zRx       = -37.5,
      times     = 10 .^LinRange(-6, -1, 50),
      ramp      = ones(10, 10),
      nfreqsperdecade = 5,
      ntimesperdecade = 10,
      glegintegorder = 5,
      lowpassfcs = [],
      provideddt = true,
      doconvramp = true,
      modelprimary = false,
      nkᵣeval = 50,
      getradialH = false,
      getazimH  = false,
      freqlow = 1e-4,
      freqhigh = 1e6,
      calcjacobian = false
  )
    @assert all(freqs .> 0.)
    @assert freqhigh > freqlow
    @assert all(diff(times) .> 0)
    thickness = zeros(nmax)
    zintfc    = zeros(nmax)
    pz        = zeros(Complex{Float64}, nmax)
    ϵᵢ      = similar(pz)
    rTE       = zeros(length(pz)-1)
    rTM       = similar(rTE)
    if freqhigh < 3/minimum(times)
       freqhigh = 3/minimum(times)
    end
    if freqlow > 0.33/maximum(times)
       freqlow = 0.33/maximum(times)
    end
    if isempty(freqs)
        freqs = 10 .^(log10(freqlow):1/nfreqsperdecade:log10(freqhigh))
    end
    interpkᵣ = 10 .^range(minimum(log10.(Filter_base))-0.5, maximum(log10.(Filter_base))+0.5, length = nkᵣeval)
    J0_kernel_h, J1_kernel_h, J0_kernel_v, J1_kernel_v = map(x->zeros(ComplexF64, length(interpkᵣ), length(freqs)), 1:4)
    J01kernelhold = zeros(ComplexF64, length(Filter_base))
    log10ω = log10.(2*pi*freqs)
    interptimes = 10 .^(minimum(log10.(times))-1:1/ntimesperdecade:maximum(log10.(times))+1)
    HFD_z       = zeros(ComplexF64, length(freqs)) # space domain fields in freq
    HFD_r       = zeros(ComplexF64, length(freqs)) # space domain fields in freq
    HFD_az       = zeros(ComplexF64, length(freqs)) # space domain fields in freq
    dBzdt     = zeros(Float64, length(times)) # time derivative of space domain fields convolved with ramp
    dBrdt     = zeros(Float64, length(times)) # time derivative of space domain fields convolved with ramp
    dBazdt     = zeros(Float64, length(times)) # time derivative of space domain fields convolved with ramp
    HFD_z_interp = zeros(ComplexF64, length(Filter_t_base))
    HFD_r_interp = zeros(ComplexF64, length(Filter_t_base))
    HFD_az_interp = zeros(ComplexF64, length(Filter_t_base))
    HTD_z_interp = zeros(Float64, length(interptimes))
    HTD_r_interp = zeros(Float64, length(interptimes))
    HTD_az_interp = zeros(Float64, length(interptimes))
    lowpassfcs = float.([lowpassfcs..., 5e6])
    quadnodes, quadweights = gausslegendre(glegintegorder)
    rxwithinloop = false
    log10Filter_base = log10.(Filter_base/rRx)
    if rTx != nothing
        if rTx>=rRx
            rxwithinloop = true
            interpkᵣ = interpkᵣ/rTx
            log10Filter_base = log10.(Filter_base/rTx)
        else
            interpkᵣ = interpkᵣ/rRx
        end
    end
    log10interpkᵣ = log10.(interpkᵣ)
    # vertical
    HFD_z_J     = calcjacobian ? zeros(ComplexF64, nmax, length(freqs)) : zeros(0, 0)
    dBzdt_J     = calcjacobian ? zeros(Float64, nmax, length(times)) : zeros(0, 0)
    HTD_z_J_interp = calcjacobian ? zeros(Float64, nmax, length(interptimes)) : zeros(0, 0)
    # azimuthal
    HFD_az_J     = calcjacobian ? zeros(ComplexF64, nmax, length(freqs)) : zeros(0, 0)
    dBazdt_J     = calcjacobian ? zeros(Float64, nmax, length(times)) : zeros(0, 0)
    HTD_az_J_interp = calcjacobian ? zeros(Float64, nmax, length(interptimes)) : zeros(0, 0)
    # radial
    HFD_r_J     = calcjacobian ? zeros(ComplexF64, nmax, length(freqs)) : zeros(0, 0)
    dBrdt_J     = calcjacobian ? zeros(Float64, nmax, length(times)) : zeros(0, 0)
    HTD_r_J_interp = calcjacobian ? zeros(Float64, nmax, length(interptimes)) : zeros(0, 0)
    #
    Jtemp = calcjacobian ? zeros(ComplexF64, nmax) : zeros(0)
    Jac = calcjacobian ? zeros(ComplexF64, nkᵣeval, nmax) : zeros(0, 0)
    useprimary = modelprimary ? one(Float64) : zero(Float64)
    HFieldDHT(thickness, pz, ϵᵢ, zintfc, rTE, rTM, zRx, zTx, rTx, rRx, freqs, times, ramp, log10ω, interptimes,
            HFD_z, HFD_r, HFD_az, HFD_z_interp, HFD_r_interp, HFD_az_interp,
            HTD_z_interp, HTD_r_interp, HTD_az_interp, dBzdt, dBrdt, dBazdt, J0_kernel_h, J1_kernel_h, J0_kernel_v, J1_kernel_v, J01kernelhold, lowpassfcs,
            quadnodes, quadweights, preallocate_ω_Hsc(interptimes, lowpassfcs)..., rxwithinloop, provideddt, doconvramp, useprimary,
            nkᵣeval, interpkᵣ, log10interpkᵣ, log10Filter_base, getradialH, getazimH, 
            calcjacobian, Jtemp, similar(Jtemp), Jac, HFD_z_J, HTD_z_J_interp, dBzdt_J,  HFD_az_J, HTD_az_J_interp, dBazdt_J, 
            HFD_r_J, HTD_r_J_interp, dBrdt_J)
end

#update geometry and dependent parameters - necessary for adjusting geometry
#e.g. inverting nuisance parameters
function update_ZR!(F::HFieldDHT, zTx, zRx, rTx, rRx)
    if F.rTx != nothing
        if F.rTx >= F.rRx
            F.interpkᵣ *= F.rTx
            F.log10interpkᵣ .+= log10(F.rTx)
            F.log10Filter_base .+= log10(F.rTx)
        else
            F.interpkᵣ *= F.rRx
            F.log10interpkᵣ .+= log10(F.rRx)
            F.log10Filter_base .+= log10(F.rRx)
        end
    else
        F.log10Filter_base .+= log10(F.rRx)
    end

    if rTx != nothing
        if rTx >= rRx
            F.interpkᵣ /= rTx
            F.log10interpkᵣ .-= log10(rTx)
            F.log10Filter_base .-= log10(rTx)
        else
            F.interpkᵣ /= rRx
            F.log10interpkᵣ .-= log10(rRx)
            F.log10Filter_base .-= log10(rRx)
        end
    else
        F.log10Filter_base .-= log10(rRx)
    end
    F.rTx = rTx
    F.rRx = rRx
    F.zTx = zTx
    F.zRx = zRx

end

function getfreqsperdec(freqvec)
    fmin, fmax = extrema(log10.(freqvec))
    ndiv = length(freqvec)-1
    freqinterval = (fmax-fmin)/ndiv
    round(1/freqinterval;sigdigits=5)
end

gettimesperdec(timevec) = getfreqsperdec(timevec) # same calculation in time

function preallocate_ω_Hsc(interptimes, lowpassfcs)
    ω   = zeros(length(Filter_t_base), length(interptimes))
    Hsc = zeros(ComplexF64, length(Filter_t_base), length(interptimes))
    for itime = 1:length(interptimes)
        t = interptimes[itime]
        ω[:,itime] = log10.(Filter_t_base/t)
        Hsc[:,itime] = ones(ComplexF64, length(Filter_t_base))
        s = 1im*10 .^ω[:,itime]
        for fc in lowpassfcs
             Hs = 1 ./( 1 .+s./(2*pi*fc))
             Hsc[:,itime] .= Hsc[:,itime].*Hs
        end
    end
    ω, 10 .^ω, Hsc
end

const μ       = 4*pi*1e-7
const ϵ₀     = 8.854e-12
const iTxLayer = 1
const iRxLayer = 1

function getpartialkz(ω, kz)
    0.5im*ω*μ/kz
end    

function getpartial_rTE(partialkz, kz, kznext)
    2*partialkz*kznext/(kz+kznext)^2
end

function getpartial_rTE_withnext(partialkznext, kz, kznext)
    -2*partialkznext*kz/(kz+kznext)^2
end    

function getpartialRstack(partial_rTE, rTE, a)
    partial_rTE*(1 - a^2)/(1 + rTE*a)^2
end 

function getpartial_bwithnext(b, dnext, partialkznext)
    b*(2im*dnext*partialkznext)
end

function stacks!(F::HField, iTxLayer::Int, nlayers::Int, ω::Float64)

    rTE              = F.rTE
    rTM              = F.rTM
    pz               = view(F.pz, 1:nlayers)
    # The first layer thicknesses is infinite
    d                = view(F.thickness, 1:nlayers)
    d[1]             = 1e60
    d[nlayers]       = 1e60
    if F.calcjacobian
        b      = view(F.b, 1:nlayers)
        J      = view(F.Jtemp, 1:nlayers)
    end    
    # Capital R is for a stack
    # Starting from the bottom up, for Rs_down
    Rlowerstack_TE, Rlowerstack_TM = zero(ComplexF64), zero(ComplexF64)
    Rlowerstack_plus_TE = zero(ComplexF64)
    @inbounds @fastmath for k = (nlayers-1):-1:iTxLayer
        if F.calcjacobian
            kz, kznext = ω*pz[k], ω*pz[k+1]
            partialkznext = getpartialkz(ω, kznext)
            S = getpartial_rTE_withnext(partialkznext, kz, kznext)
            if k == nlayers-1
                J[k+1] = S
                Rlowerstack_plus_TE = Rlowerstack_TE
                b[k+1] = exp(2im*kznext*d[k+1])
                Rlowerstack_TE = lowerstack(Rlowerstack_TE, rTE, k, b[k+1])
                continue
            end
            b[k+1] = exp(2im*kznext*d[k+1])
            partial_bnext = getpartial_bwithnext(b[k+1], d[k+1], partialkznext)
            S *= 1.0 - (Rlowerstack_TE*b[k+1])^2
            kznextplus = ω*pz[k+2]
            partial_rTE_next = getpartial_rTE(partialkznext, kznext, kznextplus)
            a = Rlowerstack_plus_TE*b[k+2]
            partialRnext = getpartialRstack(partial_rTE_next, rTE[k+1], a)
            denom = (1+rTE[k]*Rlowerstack_TE*b[k+1])^2
            S += (1.0-rTE[k]^2)*(partialRnext*b[k+1] + partial_bnext*Rlowerstack_TE) 
            S /= denom
            J[k+1] = S
            c = b[k+1]*(1.0-rTE[k]^2)/denom
            for kk = nlayers:-1:k+2
                J[kk] *= c
            end    
            Rlowerstack_plus_TE = Rlowerstack_TE
            Rlowerstack_TE = lowerstack(Rlowerstack_TE, rTE, k, b[k+1])
        else
            kznext = ω*pz[k+1]
            bnext = exp(2im*kznext*d[k+1])
            Rlowerstack_TE = lowerstack(Rlowerstack_TE, rTE, k, bnext)
        end    
    end

return Rlowerstack_TE, Rlowerstack_TM
end

@inline function lowerstack(Rlowerstack::ComplexF64, r::Array{ComplexF64, 1}, k::Int, b)
    a = Rlowerstack*b                
    Rs_d = (r[k] + a) / (1. + r[k]*a)
end

function getCurlyR(Rs_d::ComplexF64, pz::ComplexF64,
                  zR::Float64, z::SubArray{Float64, 1}, iTxLayer::Int, ω::Float64,
                  useprimary::Float64)
                  
    if (zR>=0)
        e_to_the_iwpzzr               = exp( im*ω*pz*zR)
        e_to_the_iwpz_2znext_minus_zr = exp( im*ω*pz*(2*z[iTxLayer+1] - zR))

        finRA = useprimary*e_to_the_iwpzzr + e_to_the_iwpz_2znext_minus_zr*Rs_d
        # Will need these for other components
        # finRB = e_to_the_iwpzzr + e_to_the_iwpz_2znext_minus_zr*Rs_d
        #
        finRC = useprimary*e_to_the_iwpzzr + e_to_the_iwpz_2znext_minus_zr*(-Rs_d)
        #
        finRD = useprimary*e_to_the_iwpzzr + e_to_the_iwpz_2znext_minus_zr*(-Rs_d)
    else
        e_to_the_iwpz2znext              = exp(im*ω*pz*2*z[iTxLayer+1])
        e_to_the_minus_iwpzzr            = exp(-im*ω*pz*zR)

        finRA = (useprimary + e_to_the_iwpz2znext*Rs_d) * e_to_the_minus_iwpzzr
        # Will need these for other components
        # finRB = (-1. + e_to_the_iwpz2znext*Rs_d) * e_to_the_minus_iwpzzr
        #
        finRC = (useprimary + e_to_the_iwpz2znext*(-Rs_d)) * e_to_the_minus_iwpzzr
        #
        finRD = (-useprimary + e_to_the_iwpz2znext*(-Rs_d)) * e_to_the_minus_iwpzzr
    end

    return finRA, finRD, finRC
end
# FT convention in this function means will have to conjugate for compatibility with FT = exp(-iωt)
getepsc(ρ::Float64, ω::Float64)      = Complex(ϵ₀, 1/(ρ*ω)) # corresponds to FT = exp(+iωt)
getpz(ϵᵢ::ComplexF64, kᵣ::Float64, ω::Float64) = unsafesqrt(μ*ϵᵢ - (kᵣ/ω)^2) # corresponds k² = iωμσ
ztxorignify(z, zTx)      = z - zTx
makesane(pz::Complex)    = imag(pz)  < 0.0 ? ( pz*=-1.) : pz # Sommerfeld radiation condition for fields 0 at ∞

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

function getAEM1DKernelsH!(F::HField, kᵣ::Float64, f::Float64, zz::Array{Float64, 1}, ρ::Array{Float64, 1})
    nlayers = length(ρ)
    z       = view(F.zintfc, 1:nlayers)
    zRx     = ztxorignify(F.zRx, F.zTx)
    ω       = 2. *pi*f
    ϵᵢ      = F.ϵᵢ
    pz      = F.pz
    # reflection coefficients (downward) for an intfc: pz vertical slowness
    rTE, rTM = F.rTE, F.rTM

    l = 1
    z[l]       = ztxorignify(zz[l], F.zTx)
    ϵᵢ[l]      = getepsc(ρ[l], ω)
    pz[l]      = getpz(ϵᵢ[l], kᵣ, ω)
    @inbounds @fastmath for intfc in 1:nlayers-1
        l = intfc+1
        z[l]       = ztxorignify(zz[l], F.zTx)
        ϵᵢ[l]      = getepsc(ρ[l], ω)
        pz[l]      = getpz(ϵᵢ[l], kᵣ, ω)
        rTE[intfc] = (pz[intfc] - pz[intfc+1])/(pz[intfc] + pz[intfc+1])
        # commented out as we aren't yet using TM modes
        # rTM[intfc] = (ϵᵢ[intfc]*pz[intfc+1] - ϵᵢ[intfc+1]*pz[intfc]) /
        #              (ϵᵢ[intfc]*pz[intfc+1] + ϵᵢ[intfc+1]*pz[intfc])
        F.thickness[intfc] = z[intfc+1] - z[intfc]
    end

    # TE and TM modes
    Rs_dTE, Rs_dTM   = stacks!(F, iTxLayer, nlayers, ω)
    curlyRA, curlyRD, curlyRC = getCurlyR(Rs_dTE, pz[iTxLayer], zRx, z, iTxLayer, ω, F.useprimary)
    lf_gA_TE  = 1.0
    if F.rxwithinloop
        lf_gA_TE *= loopfactor(F.rTx, kᵣ, F.rRx)
    else
        lf_gA_TE *= loopfactor(F.rTx, kᵣ)
    end
    if F.calcjacobian
        ikᵣ = findfirst(isapprox.(kᵣ, F.interpkᵣ))
        for ilayer = 2:nlayers
            curlyRAprime, = getCurlyR(F.Jtemp[ilayer], pz[iTxLayer], zRx, z, iTxLayer, ω, 0.)# cannot model primary for deriv
            F.derivmatrix[ikᵣ,ilayer] = 1/pz[iTxLayer]*curlyRAprime*lf_gA_TE*1im/(ω)*log(10)/ρ[ilayer]
        end    
    end  
    gA_TE            = μ/pz[iTxLayer]*curlyRA*lf_gA_TE
    gD_TE            = curlyRD
    gC_TE            = pz[iTxLayer]/μ*curlyRC
    # Kernels according to Loseth, without the bessel functions multiplied
    J0v              = gA_TE*1im/(ω*μ)
    J1v              = gD_TE*kᵣ^2/4pi
    J1h              = gC_TE*1im*ω*μ/4pi
    return J0v, J1v, J1h
end

loopfactor(rTx::Nothing, kᵣ::Float64)               = kᵣ^3/(4*pi)
loopfactor(rTx::Float64, kᵣ::Float64)               = kᵣ^2*besselj1(kᵣ*rTx)/(2*pi*rTx)
loopfactor(rTx::Float64, kᵣ::Float64, rRx::Float64) = kᵣ^2*besselj0(kᵣ*rRx)/(2*pi*rTx)

function getfieldFD!(F::HFieldDHT, z::Array{Float64, 1}, ρ::Array{Float64, 1})
    for (ifreq, freq) in enumerate(F.freqs)
        for (ikᵣ, kᵣ) in enumerate(F.interpkᵣ)
            if F.rxwithinloop
                F.J1_kernel_v[ikᵣ,ifreq], = getAEM1DKernelsH!(F, kᵣ, freq, z, ρ)
            else
                F.J0_kernel_v[ikᵣ,ifreq],
                F.J1_kernel_v[ikᵣ,ifreq],
                F.J1_kernel_h[ikᵣ,ifreq] = getAEM1DKernelsH!(F, kᵣ, freq, z, ρ)
            end
        end # kᵣ loop
        @views begin
            if F.rxwithinloop
                dohankeltx!(F.HFD_z, ifreq, Filter_J1, F.rTx, F.J01kernelhold, F.J1_kernel_v[:,ifreq], F.log10interpkᵣ, F.log10Filter_base)
                calcfreqdomainjacobian!(F.calcjacobian, ifreq, F.HFD_z_J, F.J01kernelhold, F.derivmatrix, F.log10interpkᵣ, F.log10Filter_base, ρ, Filter_J1, F.rTx)
            else
                # Vertical component
                dohankeltx!(F.HFD_z, ifreq, Filter_J0, F.rRx, F.J01kernelhold, F.J0_kernel_v[:,ifreq], F.log10interpkᵣ, F.log10Filter_base)
                calcfreqdomainjacobian!(F.calcjacobian, ifreq, F.HFD_z_J, F.J01kernelhold, F.derivmatrix, F.log10interpkᵣ, F.log10Filter_base, ρ, Filter_J0, F.rRx)
                # radial component
                if F.getradialH
                    dohankeltx!(F.HFD_r, ifreq, Filter_J1, F.rRx, F.J01kernelhold, F.J1_kernel_v[:,ifreq], F.log10interpkᵣ, F.log10Filter_base)
                    calcfreqdomainjacobian!(F.calcjacobian, ifreq, F.HFD_r_J, F.J01kernelhold, F.derivmatrix, F.log10interpkᵣ, F.log10Filter_base, ρ, Filter_J1, F.rRx)
                end
                # for HMD
                if F.getazimH
                    dohankeltx!(F.HFD_az, ifreq, Filter_J1, F.rRx, F.J01kernelhold, F.J1_kernel_h[:,ifreq], F.log10interpkᵣ, F.log10Filter_base)
                    calcfreqdomainjacobian!(F.calcjacobian, ifreq, F.HFD_az_J, F.J01kernelhold, F.derivmatrix, F.log10interpkᵣ, F.log10Filter_base, ρ, Filter_J1, F.rRx)
                end
            end
        end #views
    end # freq loop
end

function docomplexinterp(cmplxarray, xgiven, xwanted)
    splreal = CubicSpline(real(cmplxarray), xgiven)
    splimag = CubicSpline(imag(cmplxarray), xgiven)
    splreal.(xwanted) + 1im*splimag.(xwanted)
end    

function dohankeltx!(HFD_component, ifreq, filtercoeff, rscale, kernelhold, kernelgiven, krgiven, krwanted)
    kernelhold[:] = docomplexinterp(kernelgiven, krgiven, krwanted)
    HFD_component[ifreq] = dot(kernelhold, filtercoeff)/rscale
end

function calcfreqdomainjacobian!(doit, ifreq, HFD_component_J, kernelhold, derivmatrix, 
    log10interpkᵣ, log10Filter_base, ρ, filtercoeff, rscale)
    if doit
        for ilayer = 2:length(ρ)
            @views begin
                dohankeltx!(HFD_component_J[ilayer,:], ifreq, filtercoeff, 
                rscale, kernelhold, derivmatrix[:,ilayer], log10interpkᵣ, log10Filter_base)
            end
        end    
    end
end    

function calctimedomainjacobian(temp, HTD_component_J_interp, ilayer, itime, Filter, t, spl_component_J_real, spl_component_J_imag, l10w, w, H)
    temp[:] = -imag(conj((spl_component_J_real[ilayer].(l10w) .+ 1im*spl_component_J_imag[ilayer].(l10w)).*H)./w)*2/pi
    HTD_component_J_interp[ilayer,itime] = dot(temp, Filter)/t
end

function getfieldTD!(F::HFieldDHT, z::Array{Float64, 1}, ρ::Array{Float64, 1})
    getfieldFD!(F, z, ρ)
    spl_z_real = CubicSpline(real(F.HFD_z), F.log10ω) # TODO preallocate
    spl_z_imag = CubicSpline(imag(F.HFD_z), F.log10ω) # TODO preallocate
    if F.getradialH
        spl_r_real = CubicSpline(real(F.HFD_r), F.log10ω) # TODO preallocate
        spl_r_imag = CubicSpline(imag(F.HFD_r), F.log10ω) # TODO preallocate
    end
    if F.getazimH
        spl_az_real = CubicSpline(real(F.HFD_az), F.log10ω) # TODO preallocate
        spl_az_imag = CubicSpline(imag(F.HFD_az), F.log10ω) # TODO preallocate
    end
    if F.calcjacobian
        spl_z_J_real, spl_z_J_imag, spl_r_J_real, 
        spl_r_J_imag, spl_az_J_real, spl_az_J_imag = map(x->Vector{CubicSpline}(undef, length(ρ)), 1:6)
        for ilayer = 2:length(ρ)
            # always get vertical components
            spl_z_J_real[ilayer] = CubicSpline(real(vec(F.HFD_z_J[ilayer,:])), F.log10ω)
            spl_z_J_imag[ilayer] = CubicSpline(imag(vec(F.HFD_z_J[ilayer,:])), F.log10ω)
            # radial component
            if F.getradialH
                spl_r_J_real[ilayer] = CubicSpline(real(vec(F.HFD_r_J[ilayer,:])), F.log10ω)
                spl_r_J_imag[ilayer] = CubicSpline(imag(vec(F.HFD_r_J[ilayer,:])), F.log10ω)
            end    
            # azimuthal component
            if F.getazimH
                spl_az_J_real[ilayer] = CubicSpline(real(vec(F.HFD_az_J[ilayer,:])), F.log10ω)
                spl_az_J_imag[ilayer] = CubicSpline(imag(vec(F.HFD_az_J[ilayer,:])), F.log10ω)
            end    
        end
        temp = zeros(length(Filter_t_base))
    end        
    @views begin
        for itime = 1:length(F.interptimes)
            l10w, H = F.interplog10ω[:,itime], F.Hsc[:,itime]
            t = F.interptimes[itime]
            # Conjugate for my sign convention, apply Butterworth and inverse transform
            if F.provideddt
                # vertical
                F.HFD_z_interp[:] .= imag(conj((spl_z_real.(l10w) .+ 1im*spl_z_imag.(l10w)).*H))*2/pi
                F.HTD_z_interp[itime] = dot(F.HFD_z_interp, Filter_t_sin)/t
                if F.calcjacobian
                    for ilayer = 2:length(ρ)
                        temp[:] = imag(conj((spl_z_J_real[ilayer].(l10w) .+ 1im*spl_z_J_imag[ilayer].(l10w)).*H))*2/pi
                        F.HTD_z_J_interp[ilayer,itime] = dot(temp, Filter_t_sin)/t
                    end    
                end
                # radial
                if F.getradialH
                    F.HFD_r_interp[:] .= imag(conj((spl_r_real.(l10w) .+ 1im*spl_r_imag.(l10w)).*H))*2/pi
                    F.HTD_r_interp[itime] = dot(F.HFD_r_interp, Filter_t_sin)/t
                end
                if F.getazimH
                    F.HFD_az_interp[:] .= imag(conj((spl_az_real.(l10w) .+ 1im*spl_az_imag.(l10w)).*H))*2/pi
                    F.HTD_az_interp[itime] = dot(F.HFD_az_interp, Filter_t_sin)/t
                end
            else
                # vertical
                w = F.ω[:,itime]
                F.HFD_z_interp[:] .= -imag(conj((spl_z_real.(l10w) .+ 1im*spl_z_imag.(l10w)).*H)./w)*2/pi
                F.HTD_z_interp[itime] = dot(F.HFD_z_interp, Filter_t_cos)/t
                # radial
                if F.getradialH
                    F.HFD_r_interp[:] .= -imag(conj((spl_r_real.(l10w) .+ 1im*spl_r_imag.(l10w)).*H)./w)*2/pi
                    F.HTD_r_interp[itime] = dot(F.HFD_r_interp, Filter_t_cos)/t
                end
                if F.getazimH
                    F.HFD_az_interp[:] .= -imag(conj((spl_az_real.(l10w) .+ 1im*spl_az_imag.(l10w)).*H)./w)*2/pi
                    F.HTD_az_interp[itime] = dot(F.HFD_az_interp, Filter_t_cos)/t
                end
                if F.calcjacobian
                    for ilayer = 2:length(ρ)
                        calctimedomainjacobian(temp, HTD_z_J_interp, ilayer, itime, Filter_t_cos, t, spl_z_J_real, spl_z_J_imag, l10w, w, H)
                        calctimedomainjacobian(temp, HTD_r_J_interp, ilayer, itime, Filter_t_cos, t, spl_r_J_real, spl_r_J_imag, l10w, w, H)
                        calctimedomainjacobian(temp, HTD_az_J_interp, ilayer, itime, Filter_t_cos, t, spl_az_J_real, spl_az_J_imag, l10w, w, H)
                    end    
                end
            end
        end
    end
    if F.doconvramp
        splz = CubicSpline(F.HTD_z_interp, log10.(F.interptimes)) # TODO preallocate
        if F.getradialH
            splr = CubicSpline(F.HTD_r_interp, log10.(F.interptimes)) # TODO preallocate
        else
            splr = splz
        end
        if F.getazimH
            splaz = CubicSpline(F.HTD_az_interp, log10.(F.interptimes)) # TODO preallocate
        else
            splaz = splz
        end
        convramp!(F, splz, splr, splaz, length(ρ))
    end
end

function convramp!(F::HFieldDHT, splz::CubicSpline, splr::CubicSpline, splaz::CubicSpline, nlayers)
    fill!(F.dBzdt, 0.)
    if F.calcjacobian 
        fill!(F.dBzdt_J, 0.)
        splz_J, splr_J, splaz_J = map(x->Vector{CubicSpline}(undef, nlayers), 1:3)
        for ilayer = 2:nlayers
            splz_J[ilayer] = CubicSpline(F.HTD_z_J_interp[ilayer,:], log10.(F.interptimes)) 
            if F.getradialH
                splr_J[ilayer] = CubicSpline(F.HTD_r_J_interp[ilayer,:], log10.(F.interptimes)) 
            end
            if F.getazimH    
                splaz_J[ilayer] = CubicSpline(F.HTD_az_J_interp[ilayer,:], log10.(F.interptimes)) 
            end    
        end
        if F.getradialH
            fill!(F.dBrdt_J, 0.)
        end
        if F.getazimH
            fill!(F.dBazdt_J, 0.)
        end    
    end    
    F.getradialH && fill!(F.dBrdt, 0.)
    F.getazimH && fill!(F.dBazdt, 0.)
    for itime = 1:length(F.times)
        for iramp = 1:size(F.ramp,1)-1
            rta, rtb  = F.ramp[iramp,1], F.ramp[iramp+1,1]
            dt   = rtb - rta
            dI   = F.ramp[iramp+1,2] - F.ramp[iramp,2]
            dIdt = dI/dt

            if rta > F.times[itime]
                break
            end
            if rtb > F.times[itime] # end in this interval
                rtb = F.times[itime]
            end

            ta = F.times[itime]-rta
            tb = max(F.times[itime]-rtb, 1e-8) # rtb > rta, so make sure this is not zero...
            a, b = log10(ta), log10(tb)
            x, w = F.quadnodes, F.quadweights
            F.dBzdt[itime] += (b-a)/2*dot(getrampresponse((b-a)/2*x .+ (a+b)/2, splz), w)*dIdt
            if F.getradialH
                F.dBrdt[itime] += (b-a)/2*dot(getrampresponse((b-a)/2*x .+ (a+b)/2, splr), w)*dIdt
            end
            if F.getazimH
                F.dBazdt[itime] += (b-a)/2*dot(getrampresponse((b-a)/2*x .+ (a+b)/2, splaz), w)*dIdt
            end
            if F.calcjacobian
                for ilayer = 2:nlayers
                    F.dBzdt_J[ilayer,itime] += (b-a)/2*dot(getrampresponse((b-a)/2*x .+ (a+b)/2, splz_J[ilayer]), w)*dIdt
                    if F.getradialH
                        F.dBrdt_J[ilayer,itime] += (b-a)/2*dot(getrampresponse((b-a)/2*x .+ (a+b)/2, splr_J[ilayer]), w)*dIdt
                    end
                    if F.getazimH
                        F.dBazt_J[ilayer,itime] += (b-a)/2*dot(getrampresponse((b-a)/2*x .+ (a+b)/2, splaz_J[ilayer]), w)*dIdt
                    end    
                end    
            end
        end
    end
end

function getrampresponse(t::Array{Float64, 1}, spl::CubicSpline)
    spl.(t).*(10 .^t)*log(10)
end

end # module
