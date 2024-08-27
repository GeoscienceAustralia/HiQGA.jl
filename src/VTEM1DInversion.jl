module VTEM1DInversion
using ..AbstractOperator, ..CommonToAll
import ..AbstractOperator.get_misfit
import ..main # for McMC
import ..AbstractOperator.makeoperator
import ..AbstractOperator.getresidual # for gradientbased
import ..Model, ..Options
using ..AEM_VMD_HMD
import ..AbstractOperator.Sounding # for storing real data
import ..AbstractOperator.returnforwrite
import ..AbstractOperator.getndata
import ..AbstractOperator.plotmodelfield!
using Random, PyPlot, DelimitedFiles, LinearMaps, SparseArrays, ..GP, LinearAlgebra, Statistics

μ = AEM_VMD_HMD.μ
const pVinv = 1e12

mutable struct dBzdt<:Operator1D
    d          :: Array{Float64, 1}
    useML      :: Bool
    σ          :: Array{Float64, 1}
    F          :: AEM_VMD_HMD.HField
    z          :: Array{Float64, 1}
    nfixed     :: Int
    ρ          :: Array{Float64, 1}
    select     :: Array{Bool, 1}
    ndata      :: Int
    J          :: AbstractArray
    W          :: SparseMatrixCSC
    res        :: Vector
end

function dBzdt(;times           = [1.],
        ramp            = [0. 1],
        rTx             = 10.,  
        zTx             = -30.,
        zRx             = -30.01, # now Z down for use
        useML           = false,
        z               = [-1.],
        ρ               = [-1.],
        nfixed          = 1,
        nmax            = length(ρ)+1,
        calcjacobian    = false,
        nfreqsperdecade = 6,
        ntimesperdecade = 10,
        d               = zeros(0),
        σ               = zeros(0),
        showgates       = false,
        modelprimary    = false,
        lowpassfcs      = [],
        freqlow         = 1e-4,
        freqhigh        = 1e6, 
        )
   
    @assert size(σ)  == size(d)
    ndata, select = getndata(d)

    F = AEM_VMD_HMD.HFieldDHT(;
        times, ramp, nmax, zTx, rTx,
        rRx = 0., zRx, 
        lowpassfcs, freqlow, freqhigh,
        calcjacobian, nfreqsperdecade,
        ntimesperdecade, modelprimary
    )
    @assert length(F.thickness) >= length(z)
    # for Gauss-Newton
    res, J, W = allocateJ(F.dBzdt_J, σ, select, nfixed, length(ρ), calcjacobian)
    aem = dBzdt(d, useML, σ, F, z, nfixed, copy(ρ), select, ndata, J, W, res)
    showgates && plotwaveformgates(aem)
    aem
end

function allocateJ(FJt, σ, select, nfixed, nmodel, calcjacobian)
    if calcjacobian && !isempty(select)
        J = FJt'
        J = J[select,nfixed+1:nmodel]
    else    
        J = zeros(0)
    end
    # always return an allocated residuals and W - small price to pay I think
    # since majority of Jacobian allocations are in aem.F
    res = similar(σ[select]) 
    Wdiag = 1 ./σ[select]
    W = sparse(diagm(Wdiag))        
    return res, J, W
end

function getresidual(aem::dBzdt, log10σ::Vector{Float64}; computeJ=false)
    F = aem.F
    F.calcjacobian = computeJ
    getfield!(-log10σ, aem)
    f = F.dBzdt[aem.select]
    d = aem.d[aem.select]
    aem.res[:] = f-d
    if computeJ
        select, nfixed = aem.select, aem.nfixed
        copy!(aem.J, F.dBzdt_J'[select, nfixed+1:nfixed+length(log10σ)])
    end
    nothing    
end    

mutable struct VTEMsoundingData <: Sounding
    sounding_string :: String
    X               :: Float64
    Y               :: Float64
    Z               :: Float64
    fid             :: Float64
    linenum         :: Int
    zTx             :: Float64
    zRx             :: Float64
    rTx             :: Float64
    lowpassfcs      :: Array{Float64, 1}
    times           :: Array{Float64, 1}
    ramp            :: Array{Float64, 2}
    noise           :: Array{Float64, 1}
    data            :: Array{Float64, 1}
end

returnforwrite(s::VTEMsoundingData) = [s.X, s.Y, s.Z, s.fid, 
    s.linenum, s.zTx, s.zRx, s.rTx]

function getndata(S::VTEMsoundingData)
    getndata(S.data)[1]
end    

function VTEMsoundingData(;rRx=nothing, zRx=nothing, zTx=12.,
                            rTx=-12.,lowpassfcs=[], 
                            tx_rx_dz = -0.01, # not reading but using format Z down -ve is rx above tx
                            times=[1., 2.], ramp=[1 2; 3 4],
                            noise=[1.], data=[1.], 
                            sounding_string="sounding", X=nothing, Y=nothing, Z=nothing,
                            linenum=nothing, fid=nothing)
    @assert rTx > 0
    @assert zTx < 0 # my coordinate system z down 
    isnothing(rRx) && (rRx = 0.)
    isnothing(zRx) && (zRx = zTx + tx_rx_dz) 
    !isempty(lowpassfcs) && @assert all(lowpassfcs .> 0)
    @assert all(diff(times) .>0 )
    @assert all(diff(ramp[:,1]) .>0 )
    @assert all((noise .>0) .| isnan.(noise))
    @assert length(data) == length(noise)
    VTEMsoundingData(sounding_string, X, Y, Z, fid, linenum, zTx, zRx, rTx,
        lowpassfcs, times, ramp, noise, data)
end

function read_survey_files(;
    fname_dat="",
    fname_specs_halt="",
    frame_height = -2,
    d        = [-2, -2],
    units=1/pVinv,
    figsize = (8,4),
    fontsize = 10,
    makeqcplots = true,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.03,
    X = -1,
    Y = -1,
    Z = -1,
    fid = -1,
    linenum = -1,
    datacutoff = nothing,
    noise_scalevec = zeros(0),
    tx_rx_dz_pass_through = 0.01, # Z up GA-AEM reading convention +ve is rx above tx
    nanchar = "*")

    @assert frame_height > 0
    @assert all(d .> 0)

    @assert X > 0
    @assert Y > 0
    @assert Z > 0
    @assert linenum > 0
    @assert fid > 0
    @info "reading $fname_dat"
    soundings = readlargetextmatrix(fname_dat, startfrom, skipevery, dotillsounding)
    soundings[soundings .== nanchar] .= NaN
    easting = soundings[:,X]
    northing = soundings[:,Y]
    topo = soundings[:,Z]
    fiducial = soundings[:,fid]
    whichline = soundings[:,linenum]
    d = soundings[:,d[1]:d[2]]
    zTx = soundings[:,frame_height] # read in Z up
    zRx = -(zTx .+ tx_rx_dz_pass_through)  # my coordinate system Z down
    zTx = -zTx # my coordinate system Z down
    # check for bad z
    idxbadz = zTx .>= 0
    (sum(idxbadz) > 0 ) && @warn("kicking out $(sum(idxbadz)) bad zTx underground")
    
    @info "reading $fname_specs_halt"
    include(fname_specs_halt)
    @assert size(d, 2) == length(times)
    @assert size(d, 2) == length(σ_halt)
    σ_halt[:] .*= units
    d[:]      .*= units
    σ           = sqrt.((multnoise*d).^2 .+ (σ_halt').^2)
    if !isempty(noise_scalevec) 
        @assert length(noise_scalevec) == length(times)
        σ = σ.*noise_scalevec'
    end    
    if !isnothing(datacutoff)
        # since my dBzdt is +ve
        idxbad = d .< datacutoff
        d[idxbad] .= NaN
    end
    makeqcplots && plotsoundingdata(d, σ, times, zTx, zRx; figsize, fontsize)
    nsoundings = size(d, 1)
    s_array = Array{VTEMsoundingData, 1}(undef, nsoundings)
    fracdone = 0 
    for is in 1:nsoundings
        idxbadz[is] && continue
        l, f = Int(whichline[is]), fiducial[is]
        s_array[is] = VTEMsoundingData(;zTx=zTx[is], zRx=zRx[is], rTx, 
            times, ramp, noise=σ[is,:], data=d[is,:], lowpassfcs,
            sounding_string="sounding_$(l)_$f",
            X=easting[is], Y=northing[is], Z=topo[is], fid=f,
            linenum=l)
        fracnew = round(Int, is/nsoundings*100)
        if (fracnew-fracdone)>10
            fracdone = fracnew
            @info "read $is out of $nsoundings"
        end        
    end
    return s_array[.!idxbadz]
end

function plotsoundingdata(d, σ, times, zTx, zRx; figsize=(8,4), fontsize=1)
    f, ax = plt.subplots(2, 2, figsize=figsize, gridspec_kw=Dict("width_ratios" => [1,0.01]))
    nsoundings = size(d, 1)
    plot_d = permutedims(d)
    plot_d[plot_d .<0] .= NaN
    img = ax[1].pcolormesh(1:nsoundings, times, log10.(plot_d), shading="nearest")
    ax[1].invert_yaxis()
    cb = colorbar(img, cax=ax[1,2])
    cb.set_label("log10(dBz/dt)")
    ax[1].set_ylabel("time s")
    ax[1].tick_params(labelbottom=false)
    axx = ax[1].twiny()
    axx.semilogy(infnanmean(σ./abs.(d), 1)[:], times, "r")
    axx.semilogy(infnanmean(σ./abs.(d), 1)[:], times, "--w")
    axx.set_xlabel("avg noise fraction")
    ax[2].plot(1:nsoundings, zRx, label="zRx")
    ax[2].plot(1:nsoundings, zTx, label="zTx")
    ax[2].set_xlabel("sounding #")
    ax[2].set_ylabel("height m")
    ax[2].invert_yaxis()
    ax[2].sharex(ax[1])
    ax[2].legend()
    ax[2,2].axis("off")
    nicenup(f, fsize=fontsize)
    nothing
end

# all calling functions here for misfit, field, etc. assume model is in log10 resistivity
# SANS the top. For lower level field calculation use AEM_VMD_HMD structs

function getfield!(m::Model, aem::dBzdt)
    getfield!(m.fstar, aem)
    nothing
end

function getfield!(m::Array{Float64}, aem::dBzdt)
    copyto!(aem.ρ, aem.nfixed+1:length(aem.ρ), 10 .^m, 1:length(m))
    if (aem.ndata>0) | isempty(aem.d)
        AEM_VMD_HMD.getfieldTD!(aem.F,  aem.z, aem.ρ)
    end    
    nothing
end

function get_misfit(m::Model, opt::Options, aem::dBzdt)
    get_misfit(m.fstar, opt, aem)
end

function get_misfit(m::AbstractArray, opt::Options, aem::dBzdt)
    calcmisfit(m, opt.debug, aem)
end

function calcmisfit(m, debug, aem)
    chi2by2 = 0.0
    if !debug
        getfield!(m, aem)
        if aem.ndata>0 
            chi2by2 += getchi2by2(aem.F.dBzdt, aem.d,
                    aem.σ, aem.select, aem.useML, aem.ndata)
        end            
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
    mlfact = computeMLfactor(aem.F.dBzdt, aem.d, aem.σ, aem.select, aem.ndata)
    sqrt(mlfact)
end

# all plotting codes here assume that the model is in log10 resistivity, SANS
# the top layer resistivity. For lower level plotting use AEM_VMD_HMD structs

function plotsoundingcurve(ax, f, t; color="k", alpha=1, lw=1)
    ax.loglog(t, μ*f*pVinv, color=color, alpha=alpha, markersize=2, linewidth=lw)
end

function plotdata(ax, d, σ, t; onesigma=true)
    sigma = onesigma ? 1 : 2
    ax.errorbar(t, μ*d*pVinv, yerr = μ*sigma*pVinv*abs.(σ),
    linestyle="none", marker=".", elinewidth=0, capsize=3, color="k")
end

function plotmodelfield!(ax, iaxis, aem::dBzdt, ρ; color=nothing, alpha=1, model_lw=1, forward_lw=1)
    stepmodel(ax, iaxis, color, ρ, aem, model_lw, alpha)
    colorused = !isnothing(color) ? color : ax[iaxis].lines[end].get_color()
    getfield!(ρ, aem)
    plotsoundingcurve(ax[iaxis+1], aem.F.dBzdt, aem.F.times; color=colorused, alpha, lw=forward_lw)
end    

function initmodelfield!(aem;  onesigma=true, figsize=(8,6))
    f, ax = plt.subplots(1, 2; figsize)
    if !isempty(aem.d)
        plotdata(ax[2], aem.d, aem.σ, aem.F.times; onesigma)
    end
    ax[1].set_xlabel("log10 ρ")
    ax[1].set_ylabel("depth m")
    ax[2].set_ylabel("dBz/dt pV/Am⁴")    
    ax[2].set_xlabel("time s")
    ax[2].grid(which="both")
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
    revax && ax[1].invert_xaxis()
    nicenup(ax[1].get_figure(), fsize=12)
    ax
end 

# noisy synthetic model making
function makenoisydata!(aem, ρ;
        rseed=123, noisefrac=0.03, σ_halt=nothing, useML=false,
        onesigma=true, color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,6), revax=true,
        # σ_halt default assumed in Bfield units of pV
        units=1/pVinv)
    getfield!(ρ, aem)
    f = aem.F.dBzdt
    σ_halt = isnothing(σ_halt) ? zeros(size(f)) : units*σ_halt/μ
    σ = sqrt.((noisefrac*abs.(f)).^2 + σ_halt.^2)
    Random.seed!(rseed)
    aem.d = f + σ.*randn(size(f))
    aem.σ = σ
    aem.useML = useML
    aem.ndata, aem.select = getndata(aem.d)
    aem.res, aem.J, aem.W = allocateJ(aem.F.dBzdt_J, aem.σ, aem.select, 
        aem.nfixed, length(aem.ρ), aem.F.calcjacobian)
    plotmodelfield!(aem, ρ; onesigma, color, alpha, model_lw, forward_lw, figsize, revax)
    nothing
end

function makeoperator(sounding::VTEMsoundingData; 
            useML         = false,
            zstart        = 0.,
            zfixed        = [-1e5],
            ρfixed        = [1e12],
            extendfrac    = 1.06,
            dz            = 2.,
            ρbg           = 10, #linear ohm-m
            nlayers       = 40,
            ntimesperdecade = 10,
            modelprimary  = false,
            nfreqsperdecade = 5,
            showgeomplot  = false,
            calcjacobian  = false,
            plotfield     = false
            )
    
    zall, znall, zboundaries = setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=showgeomplot)
    z, ρ, = makezρ(zboundaries; zfixed, ρfixed)
    ρ[z.>=zstart] .= ρbg
    aem = dBzdt(;d=sounding.data/μ, σ=sounding.noise/μ, modelprimary,
        times=sounding.times, ramp=sounding.ramp, ntimesperdecade, nfreqsperdecade, lowpassfcs=sounding.lowpassfcs,
        rTx=sounding.rTx, zTx=sounding.zTx, zRx=sounding.zRx, z, ρ, calcjacobian, useML, showgates=plotfield)
    plotfield && plotmodelfield!(aem, log10.(ρ[2:end]))
    aem, zall, znall, zboundaries
end

function makeoperator(aem::dBzdt, sounding::VTEMsoundingData)
    ntimesperdecade = gettimesperdec(aem.F.interptimes)
    nfreqsperdecade = gettimesperdec(aem.F.freqs)
    modelprimary = aem.F.useprimary === 1. ? true : false
    dBzdt(;d=sounding.data/μ, σ=sounding.noise/μ, modelprimary, lowpassfcs=sounding.lowpassfcs,
        times=sounding.times, ramp=sounding.ramp, ntimesperdecade, nfreqsperdecade,
        rTx=sounding.rTx, zTx=sounding.zTx, zRx=sounding.zRx,
        z=copy(aem.z), ρ=copy(aem.ρ), 
        aem.F.calcjacobian, aem.useML, showgates=false)
end

function plotwaveformgates(aem::dBzdt; figsize=(5,5))
    figure(;figsize)
    (;ramp, times) = aem.F
    plot(ramp[:,1]*1e6, ramp[:,2], "-or")
    stem(times*1e6, ones(length(times)))
    ylabel("Amplitude")
    xlabel("time μs")
    title("Ramp and gates linear time")
end    

end
