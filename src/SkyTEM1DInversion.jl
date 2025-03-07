module SkyTEM1DInversion
import ..AbstractOperator.get_misfit, ..main
import ..AbstractOperator.Sounding
import ..AbstractOperator.makeoperator
import ..AbstractOperator.getresidual
import ..AbstractOperator.returnforwrite
import ..AbstractOperator.getndata
import ..AbstractOperator.plotmodelfield!
using ..AbstractOperator, ..AEM_VMD_HMD, Statistics, Distributed, Printf, Dates, StatsBase,
      PyPlot, LinearAlgebra, ..CommonToAll, Random, DelimitedFiles, LinearMaps, SparseArrays


import ..Model, ..Options, ..OptionsStat, ..OptionsNonstat

const μ₀ = 4*pi*1e-7
const pVinv = 1e12

mutable struct dBzdt<:Operator1D
    dlow       :: Array{Float64, 1}
    dhigh      :: Array{Float64, 1}
    useML      :: Bool
    σlow       :: Array{Float64, 1}
    σhigh      :: Array{Float64, 1}
    Flow       :: AEM_VMD_HMD.HField
    Fhigh      :: AEM_VMD_HMD.HField
    z          :: Array{Float64, 1}
    nfixed     :: Int
    ρ          :: Array{Float64, 1}
    selectlow  :: Array{Bool, 1}
    selecthigh :: Array{Bool, 1}
    ndatalow   :: Int
    ndatahigh  :: Int
    J          :: AbstractArray
    W          :: SparseMatrixCSC
    res        :: Vector
end

struct Cinv
    # assumes diagonal data cov
    oneoverσ²
end

function (a::Cinv)(x) 
    x.*a.oneoverσ²
end    

function dBzdt(;timeslow  = [.1],
                timeshigh = [.1],
                ramplow   = [.1],
                ramphigh  = [.1],
                dlow      = zeros(0),
                dhigh     = zeros(0),
                σlow      = zeros(0),
                σhigh     = zeros(0),
                rTx       = 10.,
                rRx       = 12.,
                zTxlow    = -40.,
                zTxhigh   = -40.,
                zRxlow    = -42.,
                zRxhigh   = -42.,
                nfreqsperdecade = 5,
                ntimesperdecade = 10,
                calcjacobian    = false,
                modelprimary    = false,  
                useML     = false,
                z         = [-1.],
                ρ         = [-1.],
                lowpassfcs = [],
                isRLCfilter = [],
                nfixed    = 1)
    nmax = length(ρ)+1            
    @assert size(σlow)  == size(dlow)
    @assert size(σhigh) == size(dhigh)

    Flow = AEM_VMD_HMD.HFieldDHT(;
        ntimesperdecade, nfreqsperdecade, lowpassfcs, isRLCfilter,
        times = timeslow, ramp = ramplow, nmax, zTx = zTxlow, zRx = zRxlow, 
        rRx, rTx, modelprimary, calcjacobian)

    Fhigh = AEM_VMD_HMD.HFieldDHT(;
        ntimesperdecade, nfreqsperdecade, lowpassfcs, isRLCfilter,
        times = timeshigh, ramp = ramphigh, nmax, zTx = zTxhigh, zRx = zRxhigh, 
        rRx, rTx, modelprimary, calcjacobian)

    ndatalow, selectlow    = getndata(dlow)
    ndatahigh, selecthigh  = getndata(dhigh)

    # for Gauss-Newton
    res, J, W = allocateJ(Flow, Fhigh, σlow, σhigh, selectlow, selecthigh, nfixed, length(ρ))
    dBzdt(dlow, dhigh, useML, σlow, σhigh, Flow, Fhigh, z, nfixed, copy(ρ), selectlow, selecthigh, ndatalow, ndatahigh, J, W, res)
end

function allocateJ(Flow, Fhigh, σlow, σhigh, selectlow, selecthigh, nfixed, nmodel)
    calcjacobian = Flow.calcjacobian & Fhigh.calcjacobian
    if calcjacobian && (!isempty(selectlow) || !isempty(selecthigh))
        J = [Flow.dBzdt_J'; Fhigh.dBzdt_J']; 
        J = J[[selectlow;selecthigh],nfixed+1:nmodel]
    else    
        J  = zeros(0)
    end
    # always return an allocated residuals and W - small price to pay I think
    # since majority of Jacobian allocations are in aem.F    
    Wdiag = [1 ./σlow[selectlow]; 1 ./σhigh[selecthigh]]
    res = similar(Wdiag)
    W = sparse(diagm(Wdiag))
    return res, J, W
end    

function getresidual(aem::dBzdt, log10σ::Vector{Float64}; computeJ=false)
    aem.Flow.calcjacobian = computeJ
    aem.Fhigh.calcjacobian = computeJ
    getfield!(-log10σ, aem)
    fl, fh = aem.Flow.dBzdt[aem.selectlow], aem.Fhigh.dBzdt[aem.selecthigh]
    dl, dh = aem.dlow[aem.selectlow], aem.dhigh[aem.selecthigh]
    aem.res[:] = [fl-dl; fh-dh]
    if computeJ
        selectlow, selecthigh, nfixed = aem.selectlow, aem.selecthigh, aem.nfixed
        Flow, Fhigh = aem. Flow, aem.Fhigh
        copy!(aem.J, [Flow.dBzdt_J'[selectlow,nfixed+1:nfixed+length(log10σ)]; 
                      Fhigh.dBzdt_J'[selecthigh,nfixed+1:nfixed+length(log10σ)]])
    end
    nothing    
end

mutable struct SkyTEMsoundingData <: Sounding
    sounding_string :: String
    X :: Float64
    Y :: Float64
    Z :: Float64
    fid :: Float64
    linenum :: Int
    rRx :: Float64
    zRxLM :: Float64
    zTxLM :: Float64
    zRxHM :: Float64
    zTxHM :: Float64
    rTx :: Float64
    lowpassfcs :: Array{Float64, 1}
    LM_times :: Array{Float64, 1}
    LM_ramp :: Array{Float64, 2}
    HM_times :: Array{Float64, 1}
    HM_ramp :: Array{Float64, 2}
    LM_noise :: Array{Float64, 1}
    HM_noise :: Array{Float64, 1}
    LM_data :: Array{Float64, 1}
    HM_data :: Array{Float64, 1}
    forceML
    writefields
end

returnforwrite(s::SkyTEMsoundingData) = getfield.(Ref(s), s.writefields)

function getndata(S::SkyTEMsoundingData)
    getndata(S.LM_data)[1] + getndata(S.HM_data)[1]
end   

function SkyTEMsoundingData(;rRx=-12., zRxLM=12., zTxLM=12.,
                            zRxHM=12., zTxHM=12., rTx=-12.,lowpassfcs=[-1, -2.],
                            LM_times=[1., 2.], LM_ramp=[1 2; 3 4],
                            HM_times=[1., 2.], HM_ramp=[1 2; 3 4],
                            LM_noise=[1.], HM_noise=[1.], LM_data=[1.], HM_data=[1.],
                            sounding_string="sounding", X=nothing, Y=nothing, Z=nothing,
                            linenum=nothing, fid=nothing, forceML=false, writefields=[:X, :Y, :Z, :fid, 
                            :linenum, :rRx, :zRxLM, :zTxLM, :zRxHM, :zTxHM, :rTx])
    @assert rRx > 0 && rTx > 0
    @assert zRxLM <0 && zTxLM <0
    @assert zRxHM <0 && zTxHM <0
    @assert all(lowpassfcs .> 0)
    @assert all(diff(LM_times) .>0 )
    @assert all(diff(HM_times) .>0 )
    @assert all(diff(LM_ramp[:,1]) .>0 )
    @assert all(diff(HM_ramp[:,1]) .>0 )
    @assert all((LM_noise .>0) .| isnan.(LM_noise))
    @assert all((HM_noise .>0) .| isnan.(HM_noise))
    @assert length(LM_data) == length(LM_noise)
    @assert length(HM_data) == length(HM_noise)
    SkyTEMsoundingData(sounding_string, X, Y, Z, fid, linenum, rRx, zRxLM, zTxLM, zRxHM, zTxHM, rTx,
    lowpassfcs, LM_times, LM_ramp, HM_times, HM_ramp, LM_noise, HM_noise,
    LM_data, HM_data, forceML, writefields)
end

function read_survey_files(;
    fname_dat="",
    fname_specs_halt="",
    frame_height = -2,
    frame_dz = -2,
    frame_dx = -2,
    frame_dy = -2,
    # pass through geometry takes precedence
    tx_rx_dx_pass_through = nothing, # read in Z up GA_AEM system x is fltdirn same
    tx_rx_dy_pass_through = nothing, # read in Z up GA_AEM system y has to be negated for me later down, if squared is ok
    tx_rx_dz_pass_through = nothing, # read in Z up GA_AEM system z has to be negated for me later down
    LM_Z = [-2, -2],
    HM_Z = [-2, -2],
    LM_σ = [-2, 2],
    HM_σ = [-2, 2],
    noise_scalevec = zeros(0),
    relerror = false,
    units=1/pVinv,
    figsize = (9,7),
    makeqcplots = true,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.03,
    multnoise_LM = nothing,
    multnoise_HM = nothing,
    datacutoff_LM = nothing,
    datacutoff_HM = nothing,
    X = -1,
    Y = -1,
    Z = -1,
    fid = -1,
    linenum = -1,
    nanchar = "*",
    lineslessthan = nothing,
    forceML = false, # for very low amp, in conjunction with datacutoff_LM and datacutoff_HM
    )
    @assert frame_height > 0
    @assert (frame_dz > 0) | !isnothing(tx_rx_dz_pass_through)
    @assert (frame_dx > 0) | !isnothing(tx_rx_dx_pass_through)
    @assert (frame_dy > 0) | !isnothing(tx_rx_dy_pass_through)
    @assert all(LM_Z .> 0)
    @assert all(HM_Z .> 0)
    if relerror
        @assert all(LM_σ .> 0)
        @assert all(HM_σ .> 0)
    end
    @assert X > 0
    @assert Y > 0
    @assert Z > 0
    @assert linenum > 0
    @assert fid > 0
    if !isnothing(multnoise_LM) | !isnothing(multnoise_HM)
        @assert isnothing(multnoise) "set multnoise to nothing to be consistent!"
    else
        multnoise_LM = multnoise
        multnoise_HM = multnoise
    end        
    if forceML
        @assert !isnothing(datacutoff_LM)
        @assert !isnothing(datacutoff_HM)
    end
    if !isnothing(lineslessthan)
        lineslessthan::Int
    end   
    @info "reading $fname_dat"
    soundings = readlargetextmatrix(fname_dat, startfrom, skipevery, dotillsounding)
    soundings[soundings .== nanchar] .= NaN
    easting = soundings[:,X]
    northing = soundings[:,Y]
    topo = soundings[:,Z]
    fiducial = soundings[:,fid]
    whichline = soundings[:,linenum]
    d_LM = soundings[:,LM_Z[1]:LM_Z[2]]
    d_HM = soundings[:,HM_Z[1]:HM_Z[2]]
    if relerror
        σ_LM = soundings[:,LM_σ[1]:LM_σ[2]]
        σ_HM = soundings[:,HM_σ[1]:HM_σ[2]]
    end
    zTx = soundings[:,frame_height] # read in Z up GA_AEM system
    if isnothing(tx_rx_dz_pass_through) # pass through takes precedence
        zRx = -(zTx + soundings[:,frame_dz]) # my coordinate system Z down
    else
        zRx = -(zTx .+ tx_rx_dz_pass_through)  # my coordinate system Z down
    end
    zTx = -zTx # my coordinate system Z down
    # check for bad z
    idxbadz = zTx .>= 0
    (sum(idxbadz) > 0 ) && @warn("kicking out $(sum(idxbadz)) bad zTx underground")

    if isnothing(tx_rx_dx_pass_through) # pass through takes precedence
        rRx = sqrt.(soundings[:,frame_dx].^2 + soundings[:,frame_dy].^2)
    else
        rRx = sqrt(tx_rx_dx_pass_through^2 + tx_rx_dy_pass_through^2)*ones(size(soundings, 1))
    end

    @info "reading $fname_specs_halt"
    include(fname_specs_halt)
    @assert size(d_LM, 2) == length(LM_times)
    @assert size(d_HM, 2) == length(HM_times)
    d_LM[:]     .*= units
    d_HM[:]     .*= units
    LM_noise[:] .*= units
    HM_noise[:] .*= units
    if !relerror
        @assert size(d_LM, 2) == length(LM_noise)
        @assert size(d_HM, 2) == length(HM_noise)
        σ_LM = sqrt.((multnoise_LM*d_LM).^2 .+ (LM_noise').^2)
        σ_HM = sqrt.((multnoise_HM*d_HM).^2 .+ (HM_noise').^2)
        if !isempty(noise_scalevec) 
            @assert length(noise_scalevec) == length(LM_times)+length(HM_times)
            σ_LM = σ_LM.*(noise_scalevec[1:length(LM_times)])'
            σ_HM = σ_HM.*(noise_scalevec[1:length(HM_times)])'
        end
        if !isnothing(datacutoff_LM) && !forceML
            # since my dBzdt is +ve
            idxbad = d_LM .< datacutoff_LM
            d_LM[idxbad] .= NaN
        end
        if !isnothing(datacutoff_HM) && !forceML
            # since my dBzdt is +ve
            idxbad = d_HM .< datacutoff_HM
            d_HM[idxbad] .= NaN
        end
    else
        σ_LM = sqrt.((d_LM.*σ_LM).^2 .+ (LM_noise').^2)
        σ_HM = sqrt.((d_HM.*σ_HM).^2 .+ (HM_noise').^2)
    end
    nsoundings = size(soundings, 1)
    makeqcplots && plotsoundingdata(nsoundings, LM_times, HM_times, d_LM, d_HM, 
        σ_LM, σ_HM, zRx, zTx, rRx, figsize)

    s_array = Array{SkyTEMsoundingData, 1}(undef, nsoundings)
    fracdone = 0
    countforceML = 0
    for is in 1:nsoundings
        idxbadz[is] && continue
        l, f = Int(whichline[is]), fiducial[is]
        if !isnothing(lineslessthan)
            l > lineslessthan && continue # skips high_alt and repeat lines if specified
        end  
        dlow, dhigh = vec(d_LM[is,:]), vec(d_HM[is,:])
        lowampflag = checkifdatalow(dlow, dhigh, datacutoff_LM, datacutoff_HM, forceML)
        countforceML += lowampflag ? 1 : 0
        s_array[is] = SkyTEMsoundingData(rRx=rRx[is], zRxLM=zRx[is], zTxLM=zTx[is],
            zRxHM=zRx[is], zTxHM=zTx[is], rTx=rTx, lowpassfcs=lowpassfcs,
            LM_times=LM_times, LM_ramp=LM_ramp,
            HM_times=HM_times, HM_ramp=HM_ramp,
            LM_noise=σ_LM[is,:], HM_noise=σ_HM[is,:], LM_data=dlow, HM_data=dhigh,
            sounding_string="sounding_$(l)_$f",
            X=easting[is], Y=northing[is], Z=topo[is], fid=f,
            linenum=l, forceML=lowampflag)
            fracnew = round(Int, is/nsoundings*100)
        if (fracnew-fracdone)>10
            fracdone = fracnew
            @info "read $is out of $nsoundings"
        end    
    end
    idx = [isassigned(s_array, i) for i in 1:length(s_array)]
    forceML && @info("low-amplitude forceML on $countforceML out of $(sum(idx)): $(round(100*countforceML/sum(idx)))%")
    return s_array[idx]
end

function plotsoundingdata(nsoundings, LM_times, HM_times, d_LM, d_HM, 
        LM_noise, HM_noise, zRx, zTx, rRx, figsize)
    f = figure(figsize=figsize)
    ax = Array{Any, 1}(undef, 4)
    ax[1] = subplot(2,2,1)
    plot_dLM = permutedims(d_LM)
    plot_dLM[plot_dLM .<0] .= NaN
    im1 = ax[1].pcolormesh(1:nsoundings, LM_times, log10.(plot_dLM), shading="nearest")
    ax[1].set_xlabel("sounding #")
    cbLM = colorbar(im1, ax=ax[1])
    cbLM.set_label("log10 d_LM")
    ax[1].set_ylabel("LM time s")
    axLM = ax[1].twiny()
    # axLM.semilogy(LM_noise, LM_times)
    # axLM.set_xlabel("high alt noise")
    axLM.semilogy(infnanmean(LM_noise./abs.(d_LM), 1)[:], LM_times, "r")
    axLM.semilogy(infnanmean(LM_noise./abs.(d_LM), 1)[:], LM_times, "--w")
    axLM.set_xlabel("avg LM noise fraction")
    ax[2] = subplot(2,2,3,sharex=ax[1], sharey=ax[1])
    plot_dHM = permutedims(d_HM)
    plot_dHM[plot_dHM .<0] .= NaN
    im2 = ax[2].pcolormesh(1:nsoundings, HM_times, log10.(plot_dHM), shading="nearest")
    xlabel("sounding #")
    cbHM = colorbar(im2, ax=ax[2])
    cbHM.set_label("log10 d_HM")
    ax[2].set_ylabel("HM time s")
    ax[2].invert_yaxis()
    axHM = ax[2].twiny()
    # axHM.semilogy(HM_noise, HM_times)
    # axHM.set_xlabel("high alt noise")
    axHM.semilogy(infnanmean(HM_noise./abs.(d_HM), 1)[:], HM_times, "r")
    axHM.semilogy(infnanmean(HM_noise./abs.(d_HM), 1)[:], HM_times, "--w")
    axHM.set_xlabel("avg HM noise fraction")
    ax[3] = subplot(2,2,2, sharex=ax[1])
    ax[3].plot(1:nsoundings, zRx, label="Rx")
    ax[3].plot(1:nsoundings, zTx, label="Tx")
    ax[3].legend()
    ax[3].set_xlabel("sounding #")
    ax[3].set_ylabel("height m")
    ax[3].invert_yaxis()
    ax[4] = subplot(2,2,4, sharex=ax[1])
    ax[4].plot(1:nsoundings, rRx)
    ax[4].set_xlabel("sounding #")
    ax[4].set_ylabel("rRx m")
    plt.tight_layout()
end

function checkifdatalow(d_LM, d_HM, datacutoff_LM, datacutoff_HM, forceML)
    lowflag = false
    !forceML && return lowflag
    mLM, mHM = map(xx->mean(log10.(abs.(xx))), [d_LM, d_HM])
    if (mLM < log10(datacutoff_LM)) || (mHM < log10(datacutoff_HM))
        lowflag = true
    end    
    lowflag    
end

# all calling functions here for misfit, field, etc. assume model is in log10 resistivity
# SANS the top. For lower level field calculation use AEM_VMD_HMD structs

function getfield!(m::Model, aem::dBzdt)
    getfield!(m.fstar, aem)
    nothing
end

function getfield!(m::Array{Float64}, aem::dBzdt)
    copyto!(aem.ρ, aem.nfixed+1:length(aem.ρ), 10 .^m, 1:length(m))
    ((aem.ndatalow>0) | isempty(aem.dlow)) && AEM_VMD_HMD.getfieldTD!(aem.Flow,  aem.z, aem.ρ)
    ((aem.ndatahigh>0) | isempty(aem.dhigh))  && AEM_VMD_HMD.getfieldTD!(aem.Fhigh, aem.z, aem.ρ)
    nothing
end

function get_misfit(m::Model, opt::Options, aem::dBzdt)
    get_misfit(m.fstar, opt, aem)
end

function get_misfit(m::AbstractArray, opt::Options, aem::dBzdt)
    chi2by2 = 0.0
    if !opt.debug
        getfield!(m, aem)
        aem.ndatalow>0 && (chi2by2 += getchi2by2(aem.Flow.dBzdt, aem.dlow,
                    aem.σlow, aem.selectlow, aem.useML, aem.ndatalow))
        aem.ndatahigh>0 && (chi2by2 += getchi2by2(aem.Fhigh.dBzdt, aem.dhigh,
                    aem.σhigh, aem.selecthigh, aem.useML, aem.ndatahigh))
    end
    return chi2by2
end

function getchi2by2(dBzdt, d, σ, select, useML, ndata)
    r, s, idx = dBzdt, σ, select
    r .= (r - d)./s
    if useML
        chi2by2 = 0.5*ndata*log(norm(r[idx])^2)
    else
        chi2by2 = 0.5*norm(r[idx])^2
    end
end

function computeMLfactor(dBzdt, d, σ, select, ndata)
    if ndata > 0
        r, s, idx = dBzdt, σ, select
        r .= (r - d)./s
        r[idx]'*r[idx]/ndata
    else
        NaN
    end    
end    

function computeMLfactor(aem)
    mlfact_low = computeMLfactor(aem.Flow.dBzdt, aem.dlow, aem.σlow, aem.selectlow, aem.ndatalow)
    mlfact_high = computeMLfactor(aem.Fhigh.dBzdt, aem.dhigh, aem.σhigh, aem.selecthigh, aem.ndatahigh)
    sqrt(mlfact_low), sqrt(mlfact_high)
end

function makenoisydata!(aem, ρ; 
        rseed=123, noisefrac=0.03, σ_halt_low=nothing, σ_halt_high=nothing, useML=false, showplot = true,
        onesigma=true, color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,6), revax=true,
        # σ_halt always assumed in Bfield units of pV
        units = 1/pVinv)
    getfield!(ρ, aem)
    # low moment first
    f = aem.Flow.dBzdt
    σ_halt = isnothing(σ_halt_low) ? zeros(size(f)) : units*σ_halt_low/μ₀
    σ = sqrt.((noisefrac*abs.(f)).^2 + σ_halt.^2)
    Random.seed!(rseed)
    aem.dlow = f + σ.*randn(size(f))
    aem.σlow = copy(σ)
    # high moment next
    f = aem.Fhigh.dBzdt
    σ_halt = isnothing(σ_halt_high) ? zeros(size(f)) : units*σ_halt_high/μ₀
    σ = sqrt.((noisefrac*abs.(f)).^2 + σ_halt.^2)
    aem.dhigh = f + σ.*randn(size(f))
    aem.σhigh = copy(σ)
    aem.useML = useML

    aem.ndatalow, aem.selectlow    = getndata(aem.dlow)
    aem.ndatahigh, aem.selecthigh  = getndata(aem.dhigh)

    # for Gauss-Newton
    aem.res, aem.J, aem.W = allocateJ(aem.Flow, aem.Fhigh, aem.σlow, aem.σhigh, 
                    aem.selectlow, aem.selecthigh, aem.nfixed, length(aem.ρ))
    
    if showplot
        plotwaveformgates(aem)
        plotmodelfield!(aem, ρ; onesigma, color, alpha, model_lw, forward_lw, figsize, revax)
    end
    nothing
end

function makenoisydatafile!(fname::String, aem::dBzdt, ρ::Vector{Array{Float64,1}}, xrange;
	noisefrac = 0.03, σ_halt_low=nothing, σ_halt_high=nothing, units=1/pVinv)
    d = map(zip(ρ, 1:length(ρ))) do (rho, i)
        makenoisydata!(aem, rho; 
            rseed=i, # clunky but ok
            noisefrac,  σ_halt_low,  σ_halt_high, showplot=false)
        dlow, dhigh = copy(aem.dlow), copy(aem.dhigh)     
        [i 1 xrange[i] 0 0 -aem.Flow.zTx abs(aem.Flow.zTx-aem.Flow.zRx) aem.Flow.rRx 0 dlow'*μ₀/units dhigh'*μ₀/units] 
    end
    reduce(vcat, d)
    headers = 
    """
    FID\t1
    Line\t2
    Easting\t3
    Northing\t4
    Height\t5
    frame_height\t6
    frame_dz\t7
    frame_dx\t8
    frame_dy\t9
    LM_data\t10-$(9+length(aem.dlow))
    HM_data\t$(9+length(aem.dlow)+1)-$(9+length(aem.dlow)+length(aem.dhigh)) 
    """
    f = open(fname*".hdr", "w")
    write(f, headers)
    close(f)
    writedlm(fname, d)
end

function makeoperator(sounding::SkyTEMsoundingData;
                       zfixed   = [-1e5],
                       ρfixed   = [1e12],
                       zstart = 0.0,
                       extendfrac = 1.06,
                       dz = 2.,
                       ρbg = 10,
                       nlayers = 40,
                       ntimesperdecade = 10,
                       nfreqsperdecade = 5,
                       showgeomplot = false,
                       useML = false,
                       isRLCfilter = [],
                       calcjacobian = false,
                       modelprimary = false,
                       plotfield = false)
    @assert extendfrac > 1.0
    @assert dz > 0.0
    @assert ρbg > 0.0
    @assert nlayers > 1
    nmax = nlayers+1

    zall, znall, zboundaries = setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=showgeomplot)
    z, ρ, nfixed = makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
    ρ[z.>=zstart] .= ρbg

    aem = makeoperator(sounding, ntimesperdecade, nfreqsperdecade, modelprimary, calcjacobian, useML, z, ρ, isRLCfilter)
    
    if plotfield
        plotwaveformgates(aem)
        plotmodelfield!(aem, log10.(ρ[nfixed+1:end]);)
    end

    aem, zall, znall, zboundaries
end

function makeoperator(sounding::SkyTEMsoundingData, ntimesperdecade, nfreqsperdecade, modelprimary, calcjacobian, 
        useML, z, ρ, isRLCfilter)
    dBzdt(;ntimesperdecade, nfreqsperdecade, modelprimary, calcjacobian, useML,
        timeslow = sounding.LM_times, ramplow = sounding.LM_ramp, zRxlow=sounding.zRxLM, zTxlow = sounding.zTxLM,
        timeshigh = sounding.HM_times, ramphigh = sounding.HM_ramp, zRxhigh=sounding.zRxHM, zTxhigh = sounding.zTxHM,
        rRx=sounding.rRx, rTx = sounding.rTx, z, ρ, lowpassfcs = sounding.lowpassfcs, isRLCfilter,
        dlow=sounding.LM_data/μ₀, dhigh=sounding.HM_data/μ₀, σlow=sounding.LM_noise/μ₀, σhigh=sounding.HM_noise/μ₀)
end   

function makeoperator(aem::dBzdt, sounding::SkyTEMsoundingData)
    ntimesperdecade = gettimesperdec(aem.Flow.interptimes)
    nfreqsperdecade = gettimesperdec(aem.Flow.freqs)
    @assert !(aem.useML & sounding.forceML) "useML and forceML cannot both be true"
    useML = (aem.useML | sounding.forceML) # OR logic for useML with true, true disqualified earlier
    modelprimary = aem.Flow.useprimary === 1. ? true : false
    isRLCfilter = aem.Flow.isRLCfilter[1:end-1] # the last high freq AEM_VMD_HMD_puts in, to the aem struct
    makeoperator(sounding, ntimesperdecade, nfreqsperdecade, modelprimary, aem.Flow.calcjacobian, useML, 
        copy(aem.z), copy(aem.ρ), isRLCfilter)
end

# all plotting codes here assume that the model is in log10 resistivity, SANS
# the top layer resistivity. For lower level plotting use AEM_VMD_HMD structs

function plotsoundingcurve(ax, f, t; color="k", alpha=1, lw=1)
    ax.loglog(t, μ₀*f*pVinv, color=color, alpha=alpha, markersize=2, linewidth=lw)
end

function plotdata(ax, d, σ, t; onesigma=true, dtype=:LM)
    sigma = onesigma ? 1 : 2
    label = dtype == :LM ? "low moment" : "high moment"
    color = dtype == :LM ? "g" : "m"
    ax.errorbar(t, μ₀*d*pVinv; yerr = μ₀*sigma*pVinv*abs.(σ), color,
    linestyle="none", marker=".", elinewidth=1, capsize=3, label)
end

function plotmodelfield!(ax, iaxis, aem::dBzdt, ρ; color=nothing, alpha=1, model_lw=1, forward_lw=1)
    stepmodel(ax, iaxis, color, ρ, aem, model_lw, alpha)
    getfield!(ρ, aem)
    colorused = !isnothing(color) ? color : ax[iaxis].lines[end].get_color()
    plotsoundingcurve(ax[iaxis+1], aem.Flow.dBzdt, aem.Flow.times; color=colorused, alpha, lw=forward_lw)
    plotsoundingcurve(ax[iaxis+1], aem.Fhigh.dBzdt, aem.Fhigh.times; color=colorused, alpha, lw=forward_lw)
end    

function initmodelfield!(aem;  onesigma=true, figsize=(8,6))
    f, ax = plt.subplots(1, 2; figsize)
    if !isempty(aem.dlow)
        aem.ndatalow > 0 && plotdata(ax[2], aem.dlow, aem.σlow, aem.Flow.times; onesigma, dtype=:LM)
        aem.ndatahigh > 0 && plotdata(ax[2], aem.dhigh, aem.σhigh, aem.Fhigh.times; onesigma, dtype=:HM)
    end
    ax[1].set_xlabel("log10 ρ")
    ax[1].set_ylabel("depth m")
    ax[2].set_ylabel("dBz/dt pV/Am⁴")    
    ax[2].set_xlabel("time s")
    ax
end    

function plotmodelfield!(aem::dBzdt, ρ; onesigma=true, color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,6), revax=true)
    plotmodelfield!(aem, [ρ]; onesigma, color, alpha, model_lw, forward_lw, figsize, revax) 
end  

function plotmodelfield!(aem::dBzdt, manyρ::Vector{T}; onesigma=true, 
        color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,6), revax=true) where T<:AbstractArray
    ax = initmodelfield!(aem; onesigma, figsize)
    for ρ in manyρ
        plotmodelfield!(ax, 1, aem, vec(ρ); alpha, model_lw, forward_lw, color)
    end
    ax[1].invert_yaxis()
    nicenup(ax[1].get_figure(), fsize=12)
    revax && ax[1].invert_xaxis()
    ax
end 

function plotwaveformgates(aem::dBzdt; figsize=(10,5))
    plotwaveformgates(timesLM=aem.Flow.times,  rampLM=aem.Flow.ramp,
                      timesHM=aem.Fhigh.times, rampHM=aem.Fhigh.ramp)
end

function plotwaveformgates(;timesLM=nothing, rampLM=nothing,
                            timesHM=nothing, rampHM=nothing,
                            figsize=(10,5))
    @assert timesLM != nothing
    @assert timesHM != nothing
    @assert rampLM != nothing
    @assert rampHM != nothing

    f = figure(figsize=figsize)
    s1 = subplot(121)
    plot(rampLM[:,1]*1e6, rampLM[:,2], "-or")
    stem(timesLM*1e6, ones(length(timesLM)))
    ylabel("Amplitude")
    xlabel("time μs")
    title("LM linear time")
    s2 = subplot(122, sharey=s1)
    plot(rampHM[:,1]*1e3, rampHM[:,2], "-or")
    stem(timesHM*1e3, ones(length(timesHM)))
    ylabel("Amplitude")
    xlabel("time ms")
    title("HM linear time")
    nicenup(f, fsize=12)
end

function writetabdelim(fname, opt::Options, soundings::Array{SkyTEMsoundingData, 1}; 
                        nbins=50, qp1=0.05, qp2=0.95, burninfrac=0.5,
                        zstart = 0.0,
                        extendfrac = -1,
                        dz = -1,
                        nlayers = -1,)
    @assert extendfrac > 1.0
    @assert dz > 0.0
    @assert nlayers > 1
    zall, = setupz(zstart, extendfrac, dz=dz, n=nlayers)    
    io = open(fname, "w")
    for (idx, sounding) in enumerate(soundings)
        @info "Sounding number: $idx out of $(length(soundings))"
        opt.fdataname = sounding.sounding_string*"_"
        himage, edges, CI, meanimage, meandiffimage, sdslope = make1Dhist(opt, burninfrac=burninfrac, nbins = nbins, qp1=qp1, qp2=qp2, temperaturenum=1)
        χ² = 2*assembleTat1(opt, :U, temperaturenum=1, burninfrac=0.5)
        ndata = sum(.!isnan.(sounding.LM_data)) + sum(.!isnan.(sounding.HM_data))
        χ²mean, χ²sd = mean(χ²)/ndata, std(χ²)/ndata
        towrite = [:X, :Y, :Z, :fid, :linenum, :rRx, :zRxLM, :zTxLM, :zRxHM, :zTxHM, :rTx]
        for tw in towrite # writes sounding parts
            msg = @sprintf("%e\t", getfield(sounding, tw))
            write(io, msg)
        end
        write(io, "\t")
        # zall never was missing
        for el in [vec(zall); -CI[:,3]; -CI[:,2]; -vec(meanimage); -CI[:,1]; χ²mean; χ²sd; -vec(meandiffimage); -vec(sdslope)] # concatenate everything else
            msg = @sprintf("%e\t", el)
            write(io, msg)
        end
        write(io, "\n")
        flush(io) # slower but ensures write is complete
    end
    close(io)
end    

end
