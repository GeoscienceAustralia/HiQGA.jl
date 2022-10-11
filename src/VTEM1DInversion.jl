module VTEM1DInversion
using ..AbstractOperator, ..CommonToAll
# import ..AbstractOperator.get_misfit
import ..Model, ..Options
using ..AEM_VMD_HMD
using Random, PyPlot, DelimitedFiles, LinearMaps, SparseArrays, ..GP, LinearAlgebra

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

function dBzdt(      ;times           = [1.],
                      ramp            = [0. 1],
                      rTx             = 10.,  
                      zTx             = -30.,
                      useML           = false,
                      z               = [-1.],
                      ρ               = [-1.],
                      nfixed          = 1,
                      multnoise       = 0.03,
                      nmax            = 70,
                      calcjacobian    = false,
                      nfreqsperdecade = 6,
                      ntimesperdecade = 10,
                      d               = zeros(0),
                      σ_halt          = zeros(0)
                      )
    if !@isdefined σ_halt    
        σ_halt = zeros(size(d))
    end                      
    σ = sqrt.((multnoise*abs.(d)).^2 + σ_halt.^2)
    @assert size(σ)  == size(d)
    ndata      = sum(.!isnan.(d))
    select     = .!isnan.(d)

    F = AEM_VMD_HMD.HFieldDHT(;
        times,
        ramp,
        nmax,
        zTx,
        rTx,
        rRx         = 0.,
        zRx         = zTx-0.01,
        calcjacobian,
        nfreqsperdecade,
        ntimesperdecade,
    )
    @assert length(F.thickness) >= length(z)
    # for Gauss-Newton
    res = d[select]
    if calcjacobian
        J = F.dBzdt_J'
        J = J[select,nfixed+1:length(ρ)]
        Wdiag = [1 ./σ[select]]
        res = d[select]
    else    
        res, J, Wdiag = zeros(0), zeros(0), zeros(0)
    end    
    W = sparse(diagm(Wdiag))
    dBzdt(d, useML, σ, F, z, nfixed, copy(ρ), select, ndata, J, W, res)
end

function getresidual(aem::dBzdt, log10σ::Vector{Float64}; computeJ=false)
    F = aem.F
    F.calcjacobian = computeJ
    getfield!(-log10σ, aem)
    f = F.dBzdt[aem.select]
    d = aem.d[aem.select]
    aem.res[:] = [f-d]
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
    rRx             :: Float64
    zRx             :: Float64
    zTx             :: Float64
    rTx             :: Float64
    lowpassfcs      :: Array{Float64, 1}
    times           :: Array{Float64, 1}
    ramp            :: Array{Float64, 2}
    noise           :: Array{Float64, 1}
    data            :: Array{Float64, 1}
end

function VTEMsoundingData(;rRx=nothing, zRx=nothing, zTx=12.,
                            rTx=-12.,lowpassfcs=[],
                            times=[1., 2.], ramp=[1 2; 3 4],
                            noise=[1.], data=[1.], 
                            sounding_string="sounding", X=nothing, Y=nothing, Z=nothing,
                            linenum=nothing, fid=nothing)
    @assert rTx > 0
    @assert zTx < 0 # my coordinate system z down 
    isnothing(rRx) && (rRx = 0.)
    isnothing(zRx) && (zRx = zTx-0.1) # place receiver just above tx centre
    !isempty(lowpassfc) && @assert all(lowpassfcs .> 0)
    @assert all(diff(times) .>0 )
    @assert all(diff(ramp[:,1]) .>0 )
    @assert all((noise .>0) .| isnan.(noise))
    @assert length(data) == length(noise)
    VTEMsoundingData(sounding_string, X, Y, Z, fid, linenum, rRx, zRx, zTx, rTx,
        lowpassfcs, times, ramp, noise, data)
end

function read_survey_files(;
    fname_dat="",
    fname_specs_halt="",
    frame_height = -2,
    d        = [-2, -2],
    units=1e-12,
    figsize = (8,4),
    makesounding = false,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.03,
    X = -1,
    Y = -1,
    Z = -1,
    fid = -1,
    linenum = -1,
    nanchar = "*")

    @assert frame_height > 0
    @assert all(d .> 0)

    @assert X > 0
    @assert Y > 0
    @assert Z > 0
    @assert linenum > 0
    @assert fid > 0
    @info "reading $fname_dat"
    if !isnothing(dotillsounding)
        soundings = readdlm(fname_dat)[startfrom:skipevery:dotillsounding,:]
    else
        soundings = readdlm(fname_dat)[startfrom:skipevery:end,:]
    end
    soundings[soundings .== nanchar] .= NaN
    easting = soundings[:,X]
    northing = soundings[:,Y]
    topo = soundings[:,Z]
    fiducial = soundings[:,fid]
    whichline = soundings[:,linenum]
    d = soundings[:,d[1]:d[2]]
    zTx = -soundings[:,frame_height] # my coordinate system
 
    @info "reading $fname_specs_halt"
    include(fname_specs_halt)
    @assert size(d, 2) == length(times)
    @assert size(d, 2) == length(σ_halt)
    σ_halt[:] .*= units
    d[:]      .*= units
    σ           = sqrt.((multnoise*d).^2 .+ σ_halt.^2)

    plotsoundingdata(d, σ, times, zTx, figsize)

    if makesounding
        s_array = Array{VTEMsoundingData, 1}(undef, nsoundings)
        for is in 1:nsoundings
            l, f = Int(whichline[is]), fiducial[is]
            @info "read $is out of $nsoundings"
            s_array[is] = VTEMsoundingData(;zTx=zTx[is], rTx, 
                times, ramp, noise=σ[is,:], data=d[is,:], 
                sounding_string="sounding_$(l)_$f",
                X=easting[is], Y=northing[is], Z=topo[is], fid=f,
                linenum=l)
        end
        return s_array
    end
end

function plotsoundingdata(d, σ, times, zTx, figsize)
    f, ax = plt.subplots(2, 1, figsize=figsize, sharex=true)
    nsoundings = size(soundings, 1)
    plot_d = permutedims(d)
    plot_d[plot_dLM .<0] .= NaN
    img = ax[1].pcolormesh(1:nsoundings, times, log10.(plot_d), shading="nearest")
    cbLM = colorbar(img)
    cbLM.set_label("log10(dBzdt)")
    ax[1].set_ylabel("time s")
    axx = ax[1].twiny()
    axx.semilogy(mean(σ), times)
    axx.set_xlabel("avg high alt noise")
    ax[2].plot(1:nsoundings, zTx, label="Tx")
    ax[2].set_xlabel("sounding #")
    ax[2].set_ylabel("height m")
    ax[2].legend()
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
    chi2by2 = 0.0
    if !opt.debug
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

function plotsoundingcurve(ax, f, t; color=nothing, alpha=1, lw=1)
    if isnothing(color)
        ax.loglog(t, μ*f*pVinv, alpha=alpha, markersize=2, linewidth=lw)
    else
        ax.loglog(t, μ*f*pVinv, color=color, alpha=alpha, markersize=2, linewidth=lw)
    end    
end

function plotdata(ax, d, σ, t; onesigma=true)
    sigma = onesigma ? 1 : 2
    ax.errorbar(t, μ*d, yerr = μ*sigma*abs.(σ),
    linestyle="none", marker=".", elinewidth=0, capsize=3)
end

function plotmodelfield!(ax, iaxis, aem, ρ; color=nothing, alpha=1, model_lw=1, forward_lw=1)
    nfixed = aem.nfixed
    ax[iaxis].step(ρ, aem.z[nfixed+1:end], linewidth=model_lw, alpha=alpha)
    getfield!(ρ, aem)
    plotsoundingcurve(ax[iaxis+1], aem.F.dBzdt, aem.F.times; color, alpha, lw=forward_lw)
end    

function plotmodelfield!(aem;  onesigma=true, figsize=(8,8))
    f, ax = plt.subplots(1, 2; figsize)
    if !isempty(aem.d)
        plotdata(ax[2], aem.d, aem.σ, aem.F.times; onesigma)
    end    
    ax
end    

function plotmodelfield!(aem, ρ; onesigma=true, color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,8), revax=true)
    ax = plotmodelfield!(aem; onesigma, figsize)
    plotmodelfield!(ax, 1, aem, ρ; alpha, model_lw, forward_lw, color)
    ax[1].invert_yaxis()
    revax && ax[1].invert_xaxis()
    ax
end  

function plotmodelfield!(aem, manyρ::Vector{T}; onesigma=true, 
        color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,8), revax=true) where T<:AbstractArray
    ax = plotmodelfield!(aem; onesigma, figsize)
    for ρ in manyρ
        plotmodelfield!(ax, 1, aem, ρ; alpha, model_lw, forward_lw, color)
    end
    ax[1].invert_yaxis()
    revax && ax[1].invert_xaxis()
    ax
end 

end