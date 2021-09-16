module AEM_VMD_HMD
include("DigFilters.jl")
include("TDcossin.jl")
using LinearAlgebra, SpecialFunctions, FastGaussQuadrature, DataInterpolations

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
    derivmatrix     :: Vector{Matrix}
    HFD_z_J         :: Array{ComplexF64, 2}
end

function HFieldDHT(;
      nmax      = 200,
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
    log10ω = log10.(2*pi*freqs)
    interptimes = 10 .^(minimum(log10.(times))-1:1/ntimesperdecade:maximum(log10.(times))+1)
    HFD_z       = zeros(ComplexF64, length(freqs)) # space domain fields in freq
    HFD_z_J     = zeros(ComplexF64, nmax, length(freqs))
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
    Jtemp = zeros(ComplexF64, nmax)
    Jac = map(x->zeros(ComplexF64, nmax, nkᵣeval), 1:length(freqs))
    useprimary = modelprimary ? one(Float64) : zero(Float64)
    HFieldDHT(thickness, pz, ϵᵢ, zintfc, rTE, rTM, zRx, zTx, rTx, rRx, freqs, times, ramp, log10ω, interptimes,
            HFD_z, HFD_r, HFD_az, HFD_z_interp, HFD_r_interp, HFD_az_interp,
            HTD_z_interp, HTD_r_interp, HTD_az_interp, dBzdt, dBrdt, dBazdt, J0_kernel_h, J1_kernel_h, J0_kernel_v, J1_kernel_v, lowpassfcs,
            quadnodes, quadweights, preallocate_ω_Hsc(interptimes, lowpassfcs)..., rxwithinloop, provideddt, doconvramp, useprimary,
            nkᵣeval, interpkᵣ, log10interpkᵣ, log10Filter_base, getradialH, getazimH, calcjacobian, Jtemp, Jac, HFD_z_J)
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

function getpartialRstack(partial_rTE, rTE, a)
    partial_rTE*(1 - a^2)/(1 + rTE*a)^2
end 

function getpartialRwithnext(rTE, b, Rlowerstack)
    b*(1 - rTE^2)/(1 + rTE*b*Rlowerstack)^2
end

function stacks!(F::HField, iTxLayer::Int, nlayers::Int, ω::Float64)

    rTE              = F.rTE
    rTM              = F.rTM
    pz               = view(F.pz, 1:nlayers)
    #The last and first layer thicknesses are infinite
    d                = view(F.thickness, 1:nlayers)
    d[1]             = 1e60
    d[nlayers]       = 1e60
    if F.calcjacobian
        Jtemp        = view(F.Jtemp, 1:nlayers)
    end    
    # Capital R is for a stack
    # Starting from the bottom up, for Rs_down
    Rlowerstack_TE, Rlowerstack_TM = zero(ComplexF64), zero(ComplexF64)
    @inbounds @fastmath for k = (nlayers-1):-1:iTxLayer
        if F.calcjacobian
            b = exp(2im*ω*pz[k+1]*d[k+1])
            cnext = getpartialRwithnext(rTE[k], b, Rlowerstack_TE)
            partialkz = getpartialkz(ω, im*ω*pz[k])
            partial_rTE = getpartial_rTE(partialkz, im*ω*pz[k], im*ω*pz[k+1])
            a = Rlowerstack_TE*exp(2im*ω*pz[k+1]*d[k+1]) # can write as Rlowerstack_TE*b ...
            partialRstack = getpartialRstack(partial_rTE, rTE[k], a)
            Jtemp[k+1] = partialRstack
            for j = nlayers-1:-1:k
                Jtemp[j+1] *= cnext
            end    
        end    
        Rlowerstack_TE = lowerstack(Rlowerstack_TE, pz, rTE, d, k, ω)
        #Rlowerstack_TM = lowerstack(Rlowerstack_TM, pz, rTM, d, k, ω)
    end

return Rlowerstack_TE, Rlowerstack_TM
end

@inline function lowerstack(Rlowerstack::ComplexF64, pz::SubArray{ComplexF64, 1},
                    r::Array{ComplexF64, 1}, d::SubArray{Float64, 1}, k::Int, ω::Float64)
    a = Rlowerstack*exp(2im*ω*pz[k+1]*d[k+1])
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
        ifreq = findfirst(isapprox.(ω, 2pi*F.freqs))
        ikᵣ = findfirst(isapprox.(kᵣ, F.interpkᵣ))
        for ilayer = 2:nlayers
            curlyRAprime, = getCurlyR(F.Jtemp[ilayer], pz[iTxLayer], zRx, z, iTxLayer, ω, 0.)# cannot model primary for deriv
            F.derivmatrix[ifreq][ilayer,ikᵣ] = 1/pz[iTxLayer]*curlyRAprime*lf_gA_TE*1im/(ω)*log(10)/ρ[ilayer]
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
        if F.rxwithinloop
            splreal = CubicSpline(real(F.J1_kernel_v[:,ifreq]), F.log10interpkᵣ)
            splimag = CubicSpline(imag(F.J1_kernel_v[:,ifreq]), F.log10interpkᵣ)
            J1_kernel_v = splreal.(F.log10Filter_base) + 1im*splimag.(F.log10Filter_base)
            F.HFD_z[ifreq] = dot(J1_kernel_v, Filter_J1)/F.rTx
        else
            # Vertical component
            splreal = CubicSpline(real(F.J0_kernel_v[:,ifreq]), F.log10interpkᵣ)
            splimag = CubicSpline(imag(F.J0_kernel_v[:,ifreq]), F.log10interpkᵣ)
            J0_kernel_v = splreal.(F.log10Filter_base) + 1im*splimag.(F.log10Filter_base)
            F.HFD_z[ifreq] = dot(J0_kernel_v, Filter_J0)/F.rRx
                if F.calcjacobian
                    for ilayer = 2:length(ρ)
                        splreal = CubicSpline(real(vec(F.derivmatrix[ifreq][ilayer,:])), F.log10interpkᵣ)
                        splimag = CubicSpline(imag(vec(F.derivmatrix[ifreq][ilayer,:])), F.log10interpkᵣ)
                        J0_kernel_v_prime = splreal.(F.log10Filter_base) + 1im*splimag.(F.log10Filter_base)
                        F.HFD_z_J[ilayer, ifreq] = dot(J0_kernel_v_prime, Filter_J0)/F.rRx
                    end    
                end    
            # radial component
            if F.getradialH
                splreal = CubicSpline(real(F.J1_kernel_v[:,ifreq]), F.log10interpkᵣ)
                splimag = CubicSpline(imag(F.J1_kernel_v[:,ifreq]), F.log10interpkᵣ)
                J1_kernel_v = splreal.(F.log10Filter_base) + 1im*splimag.(F.log10Filter_base)
                F.HFD_r[ifreq] = dot(J1_kernel_v, Filter_J1)/F.rRx
            end
            # for HMD
            if F.getazimH
                splreal = CubicSpline(real(F.J1_kernel_h[:,ifreq]), F.log10interpkᵣ)
                splimag = CubicSpline(imag(F.J1_kernel_h[:,ifreq]), F.log10interpkᵣ)
                J1_kernel_h = splreal.(F.log10Filter_base) + 1im*splimag.(F.log10Filter_base)
                F.HFD_az[ifreq] = dot(J1_kernel_h, Filter_J1)/F.rRx
            end
        end
    end # freq loop
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
    @views begin
        for itime = 1:length(F.interptimes)
            l10w, H = F.interplog10ω[:,itime], F.Hsc[:,itime]
            t = F.interptimes[itime]
            # Conjugate for my sign convention, apply Butterworth and inverse transform
            if F.provideddt
                # vertical
                F.HFD_z_interp[:] .= imag(conj((spl_z_real.(l10w) .+ 1im*spl_z_imag.(l10w)).*H))*2/pi
                F.HTD_z_interp[itime] = dot(F.HFD_z_interp, Filter_t_sin)/t
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
        convramp!(F, splz, splr, splaz)
    end
end

function convramp!(F::HFieldDHT, splz::CubicSpline, splr::CubicSpline, splaz::CubicSpline)
    fill!(F.dBzdt, 0.)
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
        end
    end
end

function getrampresponse(t::Array{Float64, 1}, spl::CubicSpline)
    spl.(t).*(10 .^t)*log(10)
end

end # module
