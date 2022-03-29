module SkyTEM1DInversion
import ..AbstractOperator.get_misfit, ..main, ..gradientinv
import ..AbstractOperator.Sounding
import ..AbstractOperator.makeoperator
import ..AbstractOperator.getresidual


using ..AbstractOperator, ..AEM_VMD_HMD, Statistics, Distributed, Printf, Dates,
      PyPlot, LinearAlgebra, ..CommonToAll, Random, DelimitedFiles, LinearMaps, SparseArrays, ..GP


import ..Model, ..Options, ..OptionsStat, ..OptionsNonstat

export dBzdt, plotmodelfield!, addnoise_skytem, plotmodelfield!, plotmodelfield_skytem!

const μ₀ = 4*pi*1e-7

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

function dBzdt(       Flow       :: AEM_VMD_HMD.HField,
                      Fhigh      :: AEM_VMD_HMD.HField,
                      dlow       :: Array{Float64, 1},
                      dhigh      :: Array{Float64, 1},
                      σlow       :: Array{Float64, 1},
                      σhigh      :: Array{Float64, 1};
                      useML  = false,
                      z = [-1.],
                      ρ = [-1.],
                      nfixed = 1
                    )
    @assert length(Flow.thickness) >= length(z)
    @assert length(Fhigh.thickness) >= length(z)
    @assert size(σlow)  == size(dlow)
    @assert size(σhigh) == size(dhigh)
    ndatalow   = sum(.!isnan.(dlow))
    ndatahigh  = sum(.!isnan.(dhigh))
    selectlow  = .!isnan.(dlow)
    selecthigh = .!isnan.(dhigh)
    # for Gauss-Newton
    calcjacobian = Flow.calcjacobian & Fhigh.calcjacobian
    res = [dlow[selectlow];dhigh[selecthigh]]
    if calcjacobian
        J = [Flow.dBzdt_J'; Fhigh.dBzdt_J']; 
        J = J[[selectlow;selecthigh],nfixed+1:length(ρ)]
        Wdiag = [1 ./σlow[selectlow]; 1 ./σhigh[selecthigh]]
        res = [dlow[selectlow];dhigh[selecthigh]]
    else    
        res, J, Wdiag = zeros(0), zeros(0), zeros(0)
    end    
    W = sparse(diagm(Wdiag))
    dBzdt(dlow, dhigh, useML, σlow, σhigh,
    Flow, Fhigh, z, nfixed, copy(ρ), selectlow, selecthigh, ndatalow, ndatahigh, J, W, res)
end

function getresidual(aem::dBzdt, log10σ::Vector{Float64}; computeJ=false)
    # aem.ρ[:] = 10 .^(-log10σ)
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
end

function SkyTEMsoundingData(;rRx=-12., zRxLM=12., zTxLM=12.,
                            zRxHM=12., zTxHM=12., rTx=-12.,lowpassfcs=[-1, -2.],
                            LM_times=[1., 2.], LM_ramp=[1 2; 3 4],
                            HM_times=[1., 2.], HM_ramp=[1 2; 3 4],
                            LM_noise=[1.], HM_noise=[1.], LM_data=[1.], HM_data=[1.],
                            sounding_string="sounding", X=nothing, Y=nothing, Z=nothing,
                            linenum=nothing, fid=nothing)
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
    LM_data, HM_data)
    # @show (sounding_string, rRx, zRxLM, zTxLM, zRxHM, zTxHM, rTx,
    # lowpassfcs, LM_times, LM_ramp, HM_times, HM_ramp, LM_noise, HM_noise,
    # LM_data, HM_data)
end

function read_survey_files(dfnfile::String,;
    fname_specs_halt="",
    frame_height = -2,
    frame_dz = -2,
    frame_dx = -2,
    frame_dy = -2,
    LM_Z = [-2, -2],
    HM_Z = [-2, -2],
    LM_σ = [-2, 2],
    HM_σ = [-2, 2],
    relerror = false,
    units=1e-12,
    figsize = (9,7),
    makesounding = false,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.03,
    X = -1,
    Y = -1,
    Z = -1,
    fid = -1,
    linenum = -1)

    location = findfirst(".dfn", dfnfile)[1] # get file prefix
    fname_dat = dfnfile[1:location-1]*".dat" # the name of DAT file associated with DFN
    # now read the hdr file for strings supplied in this function
    # get corresponding column numbers
    # .
    # .
    # call existing read_survey file like so

    s_array = read_survey_files(fname_dat = fname_dat,
                        fname_specs_halt = fname_specs_halt,
                        LM_Z             = LM_Z,
                        HM_Z             = HM_Z,
                        frame_height     = frame_height,
                        frame_dz         = frame_dz,
                        frame_dy         = frame_dy,
                        frame_dx         = frame_dx,
                        X                = X,
                        Y                = Y,
                        Z                = Z,
                        fid              = fid,
                        units            = units,
                        relerror         = relerror,
                        LM_σ             = LM_σ,
                        HM_σ             = HM_σ,
                        figsize          = figsize,   
                        linenum          = linenum,
                        startfrom        = startfrom,
                        skipevery        = skipevery,
                        dotillsounding   = dotillsounding,
                        multnoise        = multnoise,
                        makesounding     = makesounding)
    s_array # return sounding array
end    

function read_survey_files(;
    fname_dat="",
    fname_specs_halt="",
    frame_height = -2,
    frame_dz = -2,
    frame_dx = -2,
    frame_dy = -2,
    LM_Z = [-2, -2],
    HM_Z = [-2, -2],
    LM_σ = [-2, 2],
    HM_σ = [-2, 2],
    relerror = false,
    units=1e-12,
    figsize = (9,7),
    makesounding = false,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.03,
    X = -1,
    Y = -1,
    Z = -1,
    fid = -1,
    linenum = -1)
    @assert frame_height > 0
    @assert frame_dz > 0
    @assert frame_dx > 0
    @assert frame_dy > 0
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
    @info "reading $fname_dat"
    if dotillsounding!= nothing
        soundings = readdlm(fname_dat)[startfrom:skipevery:dotillsounding,:]
    else
        soundings = readdlm(fname_dat)[startfrom:skipevery:end,:]
    end
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
    zTx = soundings[:,frame_height]
    zRx = -(zTx + soundings[:,frame_dz])
    zTx = -zTx
    rRx = sqrt.(soundings[:,frame_dx].^2 + soundings[:,frame_dy].^2)

    @info "reading $fname_specs_halt"
    include(fname_specs_halt)
    @assert size(d_LM, 2) == length(LM_times)
    @assert size(d_HM, 2) == length(HM_times)
    if !relerror
        @assert size(d_LM, 2) == length(LM_noise)
        @assert size(d_HM, 2) == length(HM_noise)
        LM_noise[:] .*= units
        HM_noise[:] .*= units
    else
        σ_LM[:] .*= units
        σ_HM[:] .*= units
    end
    d_LM[:]     .*= units
    d_HM[:]     .*= units
    f = figure(figsize=figsize)
    ax = Array{Any, 1}(undef, 4)
    ax[1] = subplot(2,2,1)
    nsoundings = size(soundings, 1)
    plot_dLM = permutedims(d_LM)
    plot_dLM[plot_dLM .<0] .= NaN
    pcolormesh(1:nsoundings, LM_times, log10.(plot_dLM), shading="nearest")
    xlabel("sounding #")
    cbLM = colorbar()
    cbLM.set_label("log₁₀ d_LM")
    ylabel("LM time s")
    axLM = ax[1].twiny()
    axLM.semilogy(LM_noise, LM_times)
    axLM.set_xlabel("high alt noise")
    ax[2] = subplot(2,2,3,sharex=ax[1], sharey=ax[1])
    plot_dHM = permutedims(d_HM)
    plot_dHM[plot_dHM .<0] .= NaN
    pcolormesh(1:nsoundings, HM_times, log10.(plot_dHM), shading="nearest")
    xlabel("sounding #")
    cbHM = colorbar()
    cbHM.set_label("log₁₀ d_HM")
    ylabel("HM time s")
    ax[2].invert_yaxis()
    axHM = ax[2].twiny()
    axHM.semilogy(HM_noise, HM_times)
    axHM.set_xlabel("high alt noise")
    ax[3] = subplot(2,2,2, sharex=ax[1])
    plot(1:nsoundings, zRx, label="Rx")
    plot(1:nsoundings, zTx, label="Tx")
    legend()
    xlabel("sounding #")
    ylabel("height m")
    ax[3].invert_yaxis()
    ax[4] = subplot(2,2,4, sharex=ax[1])
    ax[4].plot(1:nsoundings, rRx)
    xlabel("sounding #")
    ylabel("rRx m")
    plt.tight_layout()
    if makesounding
        s_array = Array{SkyTEMsoundingData, 1}(undef, nsoundings)
        for is in 1:nsoundings
            l, f = Int(whichline[is]), fiducial[is]
            @info "read $is out of $nsoundings"
            dlow, dhigh = vec(d_LM[is,:]), vec(d_HM[is,:])
            if !relerror
                σLM = sqrt.((multnoise*dlow).^2 + LM_noise.^2)
                σHM = sqrt.((multnoise*dhigh).^2 + HM_noise.^2)
            else
                σLM, σHM = vec(σ_LM[is,:]), vec(σ_HM[is,:])
            end
            s_array[is] = SkyTEMsoundingData(rRx=rRx[is], zRxLM=zRx[is], zTxLM=zTx[is],
                zRxHM=zRx[is], zTxHM=zTx[is], rTx=rTx, lowpassfcs=lowpassfcs,
                LM_times=LM_times, LM_ramp=LM_ramp,
                HM_times=HM_times, HM_ramp=HM_ramp,
                LM_noise=σLM, HM_noise=σHM, LM_data=dlow, HM_data=dhigh,
                sounding_string="sounding_$(l)_$f",
                X=easting[is], Y=northing[is], Z=topo[is], fid=f,
                linenum=l)
        end
        return s_array
    end
end

function getfield!(m::Model, aem::dBzdt)
    getfield!(m.fstar, aem)
    nothing
end

function getfield!(m::Array{Float64}, aem::dBzdt)
    copyto!(aem.ρ, aem.nfixed+1:length(aem.ρ), 10 .^m, 1:length(m))
    aem.ndatalow>0  && AEM_VMD_HMD.getfieldTD!(aem.Flow,  aem.z, aem.ρ)
    aem.ndatahigh>0 && AEM_VMD_HMD.getfieldTD!(aem.Fhigh, aem.z, aem.ρ)
    nothing
end

function get_misfit(m::Model, opt::Options, aem::dBzdt)
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

function plotmodelfield_skytem!(Flow::AEM_VMD_HMD.HField, Fhigh::AEM_VMD_HMD.HField,
                         z::Array{Float64, 1}, ρ::Array{Float64, 1}
                        ;figsize=(10,5))
    plotwaveformgates(timesLM=Flow.times,  rampLM=Flow.ramp,
                      timesHM=Fhigh.times, rampHM=Fhigh.ramp)
    f, ax = plt.subplots(1, 2, figsize=figsize)
    ax[1].step(log10.(ρ[2:end]), z[2:end])
    ax[1].set_xlabel("log₁₀ρ")
    ax[2].set_xlabel("time s")
    ax[1].set_ylabel("depth m")
    ax[2].set_ylabel("V/Am⁴")
    AEM_VMD_HMD.getfieldTD!(Flow, z, ρ)
    AEM_VMD_HMD.getfieldTD!(Fhigh, z, ρ)
    ax[2].loglog(Flow.times,μ₀*Flow.dBzdt, label="low moment")
    ax[2].loglog(Fhigh.times,μ₀*Fhigh.dBzdt, label="high moment")
    ax[1].grid()
    ax[1].invert_yaxis()
    ax[2].grid()
    nicenup(f)
end

function addnoise_skytem(Flow::AEM_VMD_HMD.HField, Fhigh::AEM_VMD_HMD.HField,
                  z::Array{Float64, 1}, ρ::Array{Float64, 1};
                  halt_LM = nothing,
                  halt_HM = nothing,
                  noisefrac  = 0.03,
                  noisefloorlow = μ₀*1e-14,
                  noisefloorhigh = μ₀*1e-14,
                  dz = -1.,
                  extendfrac = -1.,
                  nfixed = -1,
                  figsize=(10,5),
                  rseed=42
                  )
    if halt_LM != nothing
        @assert length(halt_LM) == length(Flow.times)
    else
        halt_LM = zeros(length(Flow.times))
    end
    if halt_HM != nothing
        @assert length(halt_HM) == length(Fhigh.times)
    else
        halt_HM = zeros(length(Fhigh.times))
    end
    @assert all((nfixed, dz, extendfrac) .> 0)
    Random.seed!(rseed)
    AEM_VMD_HMD.getfieldTD!(Flow, z, ρ)
    AEM_VMD_HMD.getfieldTD!(Fhigh, z, ρ)
    dlow  = Flow.dBzdt + sqrt.((noisefrac*abs.(Flow.dBzdt)).^2 + (halt_LM/μ₀).^2).*randn(size(Flow.dBzdt))
    dhigh = Fhigh.dBzdt + sqrt.((noisefrac*abs.(Fhigh.dBzdt)).^2 + (halt_HM/μ₀).^2).*randn(size(Fhigh.dBzdt))
    dlow[abs.(dlow).<noisefloorlow] .= NaN
    dhigh[abs.(dhigh).<noisefloorhigh] .= NaN
    σlow = sqrt.((noisefrac*abs.(dlow)).^2 + (halt_LM/μ₀).^2)
    σhigh = sqrt.((noisefrac*abs.(dhigh)).^2 + (halt_HM/μ₀).^2)
    plotmodelfield!(Flow, Fhigh, z, ρ, dlow, dhigh, σlow, σhigh,
                                figsize=figsize, nfixed=nfixed,
                                dz=dz, extendfrac=extendfrac)
    # returned data is dBzdt not H if there is a μ multiplied
    return μ₀.*(dlow, dhigh, σlow, σhigh)
end

function plotmodelfield!(Flow::AEM_VMD_HMD.HField, Fhigh::AEM_VMD_HMD.HField,
                        z, ρ, dlow, dhigh, σlow, σhigh;
                        figsize=(10,5), nfixed=-1, dz=-1., extendfrac=-1.)
    # expects data and noise in units of H, i.e. B/μ
    @assert all((nfixed, dz, extendfrac) .> 0)
    f, ax = plt.subplots(1, 2, figsize=figsize)
    ax[1].step(log10.(ρ[2:end]), z[2:end])
    if dz > 0
        axn = ax[1].twinx()
        ax[1].get_shared_y_axes().join(ax[1],axn)
        axn.step(log10.(ρ[2:end]), z[2:end])
        yt = ax[1].get_yticks()[ax[1].get_yticks().>=0]
        axn.set_yticks(yt)
        axn.set_ylim(ax[1].get_ylim()[end:-1:1])
        axn.set_yticklabels(string.(Int.(round.(getn.(yt .- z[nfixed+1], dz, extendfrac)))))
    end
    ndatalow   = sum(.!isnan.(dlow))
    ndatahigh  = sum(.!isnan.(dhigh))
    if ndatalow>0
        AEM_VMD_HMD.getfieldTD!(Flow, z, ρ)
        ax[2].loglog(Flow.times,μ₀*Flow.dBzdt, label="low moment")
        ax[2].errorbar(Flow.times, μ₀*vec(dlow), yerr = μ₀*2abs.(vec(σlow)),
                            linestyle="none", marker=".", elinewidth=0, capsize=3)
    end
    if ndatahigh>0
        AEM_VMD_HMD.getfieldTD!(Fhigh, z, ρ)
        ax[2].loglog(Fhigh.times,μ₀*Fhigh.dBzdt, label="high moment")
        ax[2].errorbar(Fhigh.times, μ₀*vec(dhigh), yerr = μ₀*2abs.(vec(σhigh)),
                        linestyle="none", marker=".", elinewidth=0, capsize=3)
    end
    ax[1].grid()
    ax[2].grid()
        ax[1].step(log10.(ρ[2:end]), z[2:end])
    ax[1].set_xlabel("log₁₀ρ")
    ax[2].set_xlabel("time s")
    ax[1].set_ylabel("depth m")
    ax[2].set_ylabel("V/Am⁴")
    axn.set_ylabel("index no.")
    nicenup(f)
end

function plotmodelfield!(aem::dBzdt, Ρ::Vector{T};
                        figsize=(8,5), dz=-1., onesigma=true, onlygetMLsampled=false,
                        extendfrac=-1., fsize=10, alpha=0.1) where T<:AbstractArray
    @assert all((dz, extendfrac) .> 0)
    nfixed, z = aem.nfixed, aem.z
    if !onlygetMLsampled 
        sigma = onesigma ? 1.0 : 2.0
        f = figure(figsize=figsize)
        ax = Vector{PyPlot.PyObject}(undef, 3)
        ax[1] = subplot(121)
        ρmin, ρmax = extrema(vcat(Ρ...))
        delρ = ρmax - ρmin
        ax[1].set_xlim(ρmin-0.1delρ,ρmax+0.1delρ)
        ax[1].plot([ρmin-0.1delρ,ρmax+0.1delρ], z[nfixed+1]*[1., 1], color="b")
        ax[2] = subplot(122)
    end    
    Flow = aem.Flow
    dlow, σlow = aem.dlow, aem.σlow
    Fhigh = aem.Fhigh
    dhigh, σhigh = aem.dhigh, aem.σhigh
    if !onlygetMLsampled
        aem.ndatalow>0 && ax[2].errorbar(Flow.times, μ₀*dlow, yerr = μ₀*sigma*abs.(σlow),
                            linestyle="none", marker=".", elinewidth=0, capsize=3, label="low moment")
        aem.ndatahigh>0 && ax[2].errorbar(Fhigh.times, μ₀*dhigh, yerr = μ₀*sigma*abs.(σhigh),
                            linestyle="none", marker=".", elinewidth=0, capsize=3, label="high moment")
    else
        errorfact_low, errorfact_high = map(x->zeros(length(Ρ)), 1:2)                        
    end
    for (i, ρ) in enumerate(Ρ)
        getfield!(ρ,  aem)
        if !onlygetMLsampled
            if aem.ndatalow>0
                Flow.dBzdt[.!aem.selectlow] .= NaN
                ax[2].loglog(Flow.times,μ₀*Flow.dBzdt, "k", alpha=alpha, markersize=2)
            end
            if aem.ndatahigh>0
                Fhigh.dBzdt[.!aem.selecthigh] .= NaN
                ax[2].loglog(Fhigh.times,μ₀*Fhigh.dBzdt, "k", alpha=alpha, markersize=2)
            end
            ax[1].step(log10.(aem.ρ[2:end]), aem.z[2:end], "-k", alpha=alpha)
        else
            errorfact_low[i], errorfact_high[i] = computeMLfactor(aem)
        end    
    end
    if !onlygetMLsampled
        ax[1].grid()
        ax[1].set_ylabel("Depth m")
        ax[1].plot(xlim(), z[nfixed+1]*[1, 1], "--k")
        if dz > 0
            axn = ax[1].twinx()
            ax[1].get_shared_y_axes().join(ax[1],axn)
            yt = ax[1].get_yticks()[ax[1].get_yticks().>=z[nfixed+1]]
            axn.set_yticks(yt)
            axn.set_ylim(ax[1].get_ylim()[end:-1:1])
            axn.set_yticklabels(string.(Int.(round.(getn.(yt .- z[nfixed+1], dz, extendfrac)))))
        end
        axn.set_ylabel("Depth index", rotation=-90, labelpad=10)
        ax[1].set_xlabel("Log₁₀ρ")
        ax[1].set_title("Model")
        ax[2].set_ylabel(L"dBzdt \; V/(A.m^4)")
        ax[2].set_xlabel("Time (s)")
        ax[2].set_title("Transient response")
        ax[2].legend()
        ax[2].grid()
        ax[1].invert_xaxis()
        nicenup(f, fsize=fsize)
        nothing, nothing, nothing, nothing
    else
        mean(errorfact_low), std(errorfact_low),
        mean(errorfact_high), std(errorfact_high)
    end    
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
    ## LM operator
    Flm = AEM_VMD_HMD.HFieldDHT(
                          ntimesperdecade = ntimesperdecade,
                          nfreqsperdecade = nfreqsperdecade,
                          lowpassfcs = sounding.lowpassfcs,
                          times  = sounding.LM_times,
                          ramp   = sounding.LM_ramp,
                          nmax   = nmax,
                          zTx    = sounding.zTxLM,
                          rRx    = sounding.rRx,
                          rTx    = sounding.rTx,
                          zRx    = sounding.zRxLM,
                          modelprimary = modelprimary,
                          calcjacobian = calcjacobian
                          )
    ## HM operator
    Fhm = AEM_VMD_HMD.HFieldDHT(
                          ntimesperdecade = ntimesperdecade,
                          nfreqsperdecade = nfreqsperdecade,
                          lowpassfcs = sounding.lowpassfcs,
                          times  = sounding.HM_times,
                          ramp   = sounding.HM_ramp,
                          nmax   = nmax,
                          zTx    = sounding.zTxHM,
                          rRx    = sounding.rRx,
                          rTx    = sounding.rTx,
                          zRx    = sounding.zRxHM,
                          modelprimary = modelprimary,
                          calcjacobian = calcjacobian
                          )
    ## create operator
    dlow, dhigh, σlow, σhigh = (sounding.LM_data, sounding.HM_data, sounding.LM_noise, sounding.HM_noise)./μ₀
    aem = dBzdt(Flm, Fhm, vec(dlow), vec(dhigh), useML=useML,
                vec(σlow), vec(σhigh), z=z, ρ=ρ, nfixed=nfixed)
    if plotfield
        plotwaveformgates(aem)
        plotmodelfield!(Flm, Fhm, z, ρ, dlow, dhigh, σlow, σhigh;
                          figsize=(12,4), nfixed=nfixed, dz=dz, extendfrac=extendfrac)
    end
    aem, znall
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
    # s3 = subplot(222)
    # loglog(rampLM[:,1], rampLM[:,2], "-or")
    # stem(timesLM, ones(length(timesLM)))
    # title("LM log time")
    # s4 = subplot(224, sharex=s3, sharey=s3)
    # loglog(rampHM[:,1], rampHM[:,2], "-or")
    # stem(timesHM, ones(length(timesHM)))
    # xlabel("time s")
    # title("HM log time")
    nicenup(f, fsize=12)
end

function make_tdgp_opt(;
                    rseed = nothing,
                    znall = znall,
                    fileprefix = "sounding",
                    nmin = 2,
                    nmax = 40,
                    K = GP.OrstUhn(),
                    demean = true,
                    sdpos = 0.05,
                    sdprop = 0.05,
                    sddc = 0.008,
                    sampledc = false,
                    fbounds = [-0.5 2.5],
                    λ = [2],
                    δ = 0.1,
                    pnorm = 2,
                    save_freq = 25,
                    restart = false
                    )
    sdev_pos = [sdpos*abs(diff([extrema(znall)...])[1])]
    sdev_prop = sdprop*diff(fbounds, dims=2)[:]
    sdev_dc = sddc*diff(fbounds, dims=2)[:]
    xall = permutedims(collect(znall))
    xbounds = permutedims([extrema(znall)...])

    history_mode = "w"
	restart && (history_mode = "a")

    updatenonstat = false
    needλ²fromlog = false
    if rseed != nothing
        Random.seed!(rseed)
    end
    opt = OptionsStat(fdataname = fileprefix*"_",
                            nmin = nmin,
                            nmax = nmax,
                            xbounds = xbounds,
                            fbounds = fbounds,
                            xall = xall,
                            λ = λ,
                            δ = δ,
                            demean = demean,
                            sdev_prop = sdev_prop,
                            sdev_pos = sdev_pos,
                            sdev_dc = sdev_dc,
                            sampledc = sampledc,
                            pnorm = pnorm,
                            quasimultid = false,
                            K = K,
                            save_freq = save_freq,
                            needλ²fromlog = needλ²fromlog,
                            updatenonstat = updatenonstat,
                            dispstatstoscreen = false,
                            history_mode = history_mode
                            )
    opt
end

function makeoperatorandoptions(soundings::Array{SkyTEMsoundingData, 1};
                        rseed = nothing,
                        fileprefix = "sounding",
                        nmin = 2,
                        nmax = 40,
                        K = GP.OrstUhn(),
                        demean = true,
                        sdpos = 0.05,
                        sdprop = 0.05,
                        sddc = 0.008,
                        sampledc = false,
                        fbounds = [-0.5 2.5],
                        λ = [2],
                        δ = 0.1,
                        pnorm = 2,
                        save_freq = 25,
                        restart = false,
                        nplot = 2,
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
                        modelprimary = false
                        )
    aem, opt = [], []
    if rseed != nothing
        Random.seed!(rseed)
    end
    for idx in randperm(length(soundings))[1:nplot]
            aem, znall = makeoperator(soundings[idx],
                        zfixed = zfixed,
                        ρfixed = ρfixed,
                        zstart = zstart,
                        extendfrac = extendfrac,
                        dz = dz,
                        ρbg = ρbg,
                        nlayers = nlayers,
                        ntimesperdecade = ntimesperdecade,
                        nfreqsperdecade = nfreqsperdecade,
                        useML = useML,
                        showgeomplot = showgeomplot,
                        modelprimary = modelprimary,
                        plotfield = true)

                    opt = make_tdgp_opt(znall = znall,
                        fileprefix = soundings[idx].sounding_string,
                        nmin = nmin,
                        nmax = nmax,
                        K = K,
                        demean = demean,
                        sdpos = sdpos,
                        sdprop = sdprop,
                        sddc = sddc,
                        sampledc = sampledc,
                        fbounds = fbounds,
                        save_freq = save_freq,
                        λ = λ,
                        δ = δ,
                        restart = restart
                        )
    end
    # returns generically all the required elements in operator
    # for plotting except sounding_string
    opt
end

function summarypost(soundings::Array{SkyTEMsoundingData, 1}, opt::Options;
            qp1=0.05,
            qp2=0.95,
            burninfrac=0.5,
            zstart = 0.0,
            extendfrac = -1,
            dz = -1,
            nlayers = -1,
            useML=false)

    @assert extendfrac > 1.0
    @assert dz > 0.0
    @assert nlayers > 1
    zall, = setupz(zstart, extendfrac, dz=dz, n=nlayers)

    linename = "_line_$(soundings[1].linenum)_summary.txt"
    fnames = ["rho_low", "rho_mid", "rho_hi", "rho_avg",
              "ddz_mean", "ddz_sdev", "phid_mean", "phid_sdev"].*linename
    if isfile(fnames[1])
        @warn fnames[1]*" exists, reading stored values"
        pl, pm, ph, ρmean,
        vdmean, vddev, χ²mean, χ²sd, = map(x->readdlm(x), fnames)
    else
        # this is a dummy operator for plotting
        aem, = makeoperator(soundings[1])
        # now get the posterior marginals
        opt.xall[:] .= zall
        pl, pm, ph, ρmean, vdmean, vddev = map(x->zeros(length(zall), length(soundings)), 1:6)
        χ²mean, χ²sd = zeros(length(soundings)), zeros(length(soundings))
        for idx = 1:length(soundings)
            @info "$idx out of $(length(soundings))\n"
            opt.fdataname = soundings[idx].sounding_string*"_"
            opt.xall[:] .= zall
            pl[:,idx], pm[:,idx], ph[:,idx], ρmean[:,idx],
            vdmean[:,idx], vddev[:,idx] = CommonToAll.plot_posterior(aem, opt, burninfrac=burninfrac,
                                                    qp1=qp1, qp2=qp2,
                                                    doplot=false)
            χ² = 2*CommonToAll.assembleTat1(opt, :U, temperaturenum=1, burninfrac=burninfrac)
            ndata = sum(.!isnan.(soundings[idx].LM_data)) +
                    sum(.!isnan.(soundings[idx].HM_data))
            χ²mean[idx] = mean(χ²)/ndata
            χ²sd[idx]   = std(χ²)/ndata
            useML && (χ²mean[idx] -= log(ndata))
        end
        # write in grid format
        for (fname, vals) in Dict(zip(fnames, [pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd]))
            writedlm(fname, vals)
        end
        # write in x, y, z, rho format
        for (i, d) in enumerate([pl, pm, ph, ρmean])
            xyzrho = makearray(soundings, d, zall)
            writedlm(fnames[i][1:end-4]*"_xyzrho.txt", xyzrho)
        end
    end
    pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd, zall
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
    end
    close(io)
end    

function plotindividualsoundings(soundings::Array{SkyTEMsoundingData, 1}, opt::Options;
    burninfrac=0.5,
    nbins=100,
    figsize  = (12,6),
    zfixed   = [-1e5],
    ρfixed   = [1e12],
    zstart = 0.0,
    extendfrac = 1.06,
    dz = 2.,
    ρbg = 10,
    omittemp = false,
    nlayers = 40,
    ntimesperdecade = 10,
    nfreqsperdecade = 5,
    computeforwards = false,
    nforwards = 100,
    rseed = 11,
    modelprimary = false,
    idxcompute = [1],
    onlygetMLsampled = false)
    
    if onlygetMLsampled
        computeforwards = true
        errorfac_low, errorfacσ_low, errorfac_high, errorfacσ_high = map(x->zeros(length(soundings)), 1:4)
    end
    zall, = setupz(zstart, extendfrac, dz=dz, n=nlayers)
    opt.xall[:] .= zall
    for idx = 1:length(soundings)
        if in(idx, idxcompute)
            @info "Sounding number: $idx"
            opt.fdataname = soundings[idx].sounding_string*"_"
            aem, znall = makeoperator(soundings[idx],
                zfixed = zfixed,
                ρfixed = ρfixed,
                zstart = zstart,
                extendfrac = extendfrac,
                dz = dz,
                ρbg = ρbg,
                nlayers = nlayers,
                ntimesperdecade = ntimesperdecade,
                nfreqsperdecade = nfreqsperdecade,
                modelprimary = modelprimary)
            if !onlygetMLsampled
                getchi2forall(opt, alpha=0.8, omittemp=omittemp)
                CommonToAll.getstats(opt)
                plot_posterior(aem, opt, burninfrac=burninfrac, nbins=nbins, figsize=figsize)
                ax = gcf().axes
                ax[1].invert_xaxis()
            end    
            if computeforwards
                M = assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)
                Random.seed!(rseed)
                if onlygetMLsampled
                    errorfac_low[idx], errorfacσ_low[idx], 
                    errorfac_high[idx], errorfacσ_high[idx] = plotmodelfield!(aem, M[randperm(length(M))[1:nforwards]],
                        dz=dz, extendfrac=extendfrac, onesigma=false, alpha=0.2, onlygetMLsampled=true)
                else
                    plotmodelfield!(aem, M[randperm(length(M))[1:nforwards]],
                    dz=dz, extendfrac=extendfrac, onesigma=false, alpha=0.2)
                end            
            end
        end
    end
    if onlygetMLsampled
        errorfac_low, errorfacσ_low, errorfac_high, errorfacσ_high
    end    
end

function plotsummarygrids1(meangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname; qp1=0.05, qp2=0.95,
                        figsize=(10,10), fontsize=12, cmap="viridis", vmin=-2, vmax=0.5, Eislast=true, Nislast=true,
                        topowidth=2, idx=nothing, omitconvergence=false, useML=false, preferEright=false, preferNright=false,
                        saveplot=false, yl=nothing, botadjust=[0.125, 0.05, 0.75, 0.01], dpi=300, showplot=true)
    f = figure(figsize=figsize)
    dr = diff(gridx)[1]
    f.suptitle(lname*" Δx=$dr m, Fids: $(length(R))", fontsize=fontsize)
    nrows = omitconvergence ? 4 : 5
    icol = 1
    s = Array{Any, 1}(undef, nrows)
    if !omitconvergence
        s[icol] = subplot(nrows, 1, icol)
        if useML
            plot(R, exp.(χ²mean))
            semilogy(R, ones(length(R)), "--k")
            ylabel("variance factor")
            title("Max likelihood variance adjustment")
        else
            plot(R, χ²mean)
            plot(R, ones(length(R)), "--k")
            fill_between(R, vec(χ²mean-χ²sd), vec(χ²mean+χ²sd), alpha=0.5)
            ylabel(L"ϕ_d")
            title("Data misfit")
        end
        icol += 1
    end
    s[icol] = omitconvergence ? subplot(nrows, 1, icol) : subplot(nrows, 1, icol, sharex=s[icol-1])
    imshow(plgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    title("Percentile $(round(Int, 100*qp1)) conductivity")
    ylabel("Height m")
    icol += 1
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1], sharey=s[icol-1])
    imshow(pmgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("Percentile 50 conductivity")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    ylabel("Height m")
    icol += 1
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1], sharey=s[icol-1])
    imshow(meangrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("Mean conductivity")
    ylabel("Height m")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    icol +=1
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1], sharey=s[icol-1])
    imlast = imshow(phgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    xlabel("Line distance m")
    title("Percentile $(round(Int, 100*qp2)) conductivity")
    ylabel("Height m")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    xlim(extrema(gridx))
    map(x->x.tick_params(labelbottom=false), s[1:end-1])
    map(x->x.grid(), s)
    nicenup(f, fsize=fontsize)
    isa(yl, Nothing) || s[end].set_ylim(yl...)
    plotNEWSlabels(Eislast, Nislast, gridx, gridz, s)
    f.subplots_adjust(bottom=0.125)
    cbar_ax = f.add_axes(botadjust)
    cb = f.colorbar(imlast, cax=cbar_ax, orientation="horizontal",)
    cb.ax.set_xlabel("Log₁₀ S/m", fontsize=fontsize)
    cb.ax.tick_params(labelsize=fontsize)
    (preferNright && !Nislast) && s[end].invert_xaxis()
    (preferEright && !Eislast) && s[end].invert_xaxis()
    saveplot && savefig(lname*".png", dpi=dpi)
    showplot || close(f)
end

function plotsummarygrids2(σmeangrid, ∇zmeangrid, ∇zsdgrid, cigrid, gridx, gridz, topofine, lname;
        qp1=0.05, qp2=0.95, Eislast=true, Nislast=true,
        figsize=(10,10), fontsize=12, cmap="viridis", vmin=-2, vmax=0.5, topowidth=2)
    f = figure(figsize=figsize)
    f.suptitle(lname, fontsize=fontsize)
    s1 = subplot(411)
    imshow(σmeangrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    title("Mean conductivity")
    ylabel("Height m")
    colorbar()
    s2 = subplot(412, sharex=s1)
    imshow(abs.(cigrid), cmap=cmap, aspect="auto",
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("CI: $(round(Int, 100*(qp2-qp1))) conductivity")
    ylabel("Height m")
    colorbar()
    s3 = subplot(413, sharex=s2, sharey=s2)
    imshow(∇zmeangrid, cmap=cmap, aspect="auto",
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("Mean conductivity vertical derivative")
    ylabel("Height m")
    colorbar()
    s4 = subplot(414, sharex=s2, sharey=s2)
    imlast = imshow(∇zsdgrid, cmap=cmap, aspect="auto",
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    xlabel("Line distance m")
    title("Std dev of conductivity vertical derivative")
    ylabel("Height m")
    colorbar()
    xlim(extrema(gridx))
    nicenup(f, fsize=fontsize)
end   

function splitlinesummaryimages(soundings::Array{SkyTEMsoundingData, 1}, opt::Options;
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        zstart = 0.0,
                        extendfrac = -1,
                        dz = -1,
                        dr = 10,
                        nlayers = -1,
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap="viridis",
                        figsize=(6,10),
                        topowidth=2,
                        idx = nothing,
                        showderivs = false,
                        omitconvergence = false,
                        useML = false,
                        preferEright = false,
                        preferNright = false,
                        saveplot = true,
                        yl = nothing,
                        showplot = true,
                        botadjust=[0.125, 0.05, 0.75, 0.01],
                        dpi = 300)
    
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)                   
    for i in 1:nlines
        a = linestartidx[i]
        b = i != nlines ?  linestartidx[i+1]-1 : length(soundings)
        summaryimages(soundings[a:b], opt, qp1=qp1, qp2=qp2, burninfrac=burninfrac, zstart=zstart,
            extendfrac=extendfrac, dz=dz, dr=dr, nlayers=nlayers, fontsize=fontsize, vmin=vmin, 
            vmax = vmax, cmap=cmap, figsize=figsize, topowidth=topowidth, idx=idx, showderivs=showderivs,
            omitconvergence=omitconvergence, useML=useML, preferEright=preferEright, showplot=showplot,
            preferNright=preferNright, saveplot=saveplot, yl=yl, botadjust=botadjust, dpi=dpi)
    end
    nothing    
end

function summaryimages(soundings::Array{SkyTEMsoundingData, 1}, opt::Options;
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        zstart = 0.0,
                        extendfrac = -1,
                        dz = -1,
                        dr = 10,
                        nlayers = -1,
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap="viridis",
                        figsize=(6,10),
                        topowidth=2,
                        idx = nothing,
                        showderivs = false,
                        omitconvergence = false,
                        useML = false,
                        preferEright = false,
                        preferNright = false,
                        saveplot = false,
                        yl = nothing,
                        showplot = true,
                        botadjust=[0.125, 0.05, 0.75, 0.01],
                        dpi = 300)
    @assert !(preferNright && preferEright) # can't prefer both labels to the right
    pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd, zall = summarypost(soundings, opt,
                                                                    zstart=zstart,
                                                                    qp1=qp1,
                                                                    qp2=qp2,
                                                                    dz=dz,
                                                                    extendfrac=extendfrac,
                                                                    nlayers=nlayers,
                                                                    burninfrac=burninfrac, useML=useML)

    phgrid, plgrid, pmgrid, σmeangrid, ∇zmeangrid,
    ∇zsdgrid, gridx, gridz, topofine, R, Z = makesummarygrid(soundings, pl, pm, ph, ρmean,
                                                            vdmean, vddev, zall, dz, dr=dr)

    lname = "Line $(soundings[1].linenum)"
    Eislast, Nislast = whichislast(soundings)
    plotsummarygrids1(σmeangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname, qp1=qp1, qp2=qp2,
                        figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax, Eislast=Eislast,
                        Nislast=Nislast, topowidth=topowidth, idx=idx, omitconvergence=omitconvergence, useML=useML,
                        preferEright=preferEright, preferNright=preferNright, saveplot=saveplot, showplot=showplot, dpi=dpi,
                        yl=yl, botadjust=botadjust)                  
    if showderivs
        cigrid = phgrid - plgrid
        plotsummarygrids2(σmeangrid, ∇zmeangrid, ∇zsdgrid, cigrid, gridx, gridz, topofine, lname, qp1=qp1, qp2=qp2,
                        figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax, topowidth=topowidth,
                        Eislast=Eislast, Nislast=Nislast)
    end
end

# for deterministic inversions, read in
function readingrid(soundings, zall)
    nstepsmax = 0 
    for s in soundings
        fname = s.sounding_string*"_gradientinv.dat"
        nstepsmax = max(size(readdlm(fname), 1), nstepsmax)
    end
    ϕd = zeros(length(soundings))
    σgrid = zeros(length(zall), length(soundings))
    for (i, s) in enumerate(soundings)
        fname = s.sounding_string*"_gradientinv.dat"
        A = readdlm(fname)
        ϕd[i] = A[end,2]
        σgrid[:,i] = vec(A[end,3:end])
    end    
    ϕd, σgrid
end

# plot the convergence and the result
function splitlineconvandlast(soundings, delr, delz; 
        zstart=-1, extendfrac=-1, dz=-1, nlayers=-1, cmapσ="jet", vmin=-2.5, vmax=0.5, fontsize=12,
        figsize=(20,5),
        topowidth=1,
        preferEright = false,
        preferNright = false,
        saveplot = true,
        showplot = true,
        dpi=400)
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)                   
    for i in 1:nlines
        a = linestartidx[i]
        b = i != nlines ?  linestartidx[i+1]-1 : length(soundings)
        plotconvandlast(soundings[a:b], delr, delz, 
            zstart=zstart, extendfrac=extendfrac, dz=dz, nlayers=nlayers, 
            cmapσ=cmapσ, vmin=vmin, vmax=vmax, fontsize=fontsize,
            figsize=figsize, topowidth=topowidth, preferEright=preferEright,
            preferNright=preferNright, saveplot=saveplot, showplot=showplot, dpi=dpi)
    end
    nothing
end

function plotconvandlast(soundings, delr, delz; 
        zstart=-1, extendfrac=-1, dz=-1, nlayers=-1, cmapσ="jet", vmin=-2.5, vmax=0.5, fontsize=12,
        figsize=(20,5),
        topowidth=1,
        preferEright = false,
        preferNright = false,
        saveplot = true, 
        showplot = true,
        dpi = 400)
    @assert zstart > -1
    @assert extendfrac >-1
    @assert dz>-1
    @assert nlayers>-1
    Eislast, Nislast = whichislast(soundings)
    zall, = setupz(zstart, extendfrac, dz=dz, n=nlayers)
    ϕd, σ = readingrid(soundings, zall)
    img, gridr, gridz, topofine, R = makegrid(σ, soundings, zall=zall, dz=delz, dr=delr)
    f, ax = plt.subplots(2,1, sharex=true, figsize=figsize)
    lname = "Line_$(soundings[1].linenum)"
    f.suptitle(lname*" Δx=$delr m, Fids: $(length(R))", fontsize=fontsize)
    ax[1].plot(R, ϕd)
    # ax[1].plot(R, ones(length(ϕd)), "--k")
    ax[1].set_ylabel(L"\phi_d")
    ax[1].set_yscale("log")
    imlast = ax[2].imshow(img, extent=[gridr[1], gridr[end], gridz[end], gridz[1]], cmap=cmapσ, aspect="auto", vmin=vmin, vmax=vmax)
    ax[2].plot(gridr, topofine, linewidth=topowidth, "-k")
    ax[2].set_xlim(extrema(gridr)...)
    ax[2].set_ylabel("mAHD")
    ax[2].set_xlabel("Distance m")
    nicenup(f, fsize=fontsize)
    plotNEWSlabels(Eislast, Nislast, gridr, gridz, [ax[2]])
    f.subplots_adjust(bottom=0.25)
    cbar_ax = f.add_axes([0.3, 0.1, 0.4, 0.02])
    cb = f.colorbar(imlast, cax=cbar_ax, orientation="horizontal")
    cb.ax.set_xlabel("Log₁₀ S/m")
    (preferNright && !Nislast) && ax[end].invert_xaxis()
    (preferEright && !Eislast) && ax[end].invert_xaxis()
    saveplot && savefig(lname*".png", dpi=dpi)
    showplot || close(f)
end    

# plot multiple grids with supplied labels
function plotgrids()
    f, ax = plt.subplots(size(grids, 1), 1, figsize=figsize,
                        sharex=true, sharey=true, squeeze=false)
    for (i, g) in enumerate(grids)
        ax[i].imshow(g, cmap=cmap, aspect="auto", alpha=α, vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
        colorbar(img, ax=ax[i])
        titles == nothing || ax[i].set_title(titles[i])
        idx  == nothing || plotprofile.(ax[i], idx)
        elev == nothing || plot(gridx, elev, "-k")
        ylabels == nothing || ax[i].set_ylabel(ylabels[i])
    end
    ax[end].set_xlabel("Easting m")
    transD_GP.CommonToAll.nicenup(f, fsize=fsize)
end

# no nuisances e.g., SkyTEM, transD
function loopacrosssoundings(soundings::Array{S, 1}, opt_in::Options;
                            nsequentialiters   =-1,
                            nparallelsoundings =-1,
                            zfixed             = [-1e5],
                            ρfixed             = [1e12],
                            useML              = false,
                            zstart             = 0.0,
                            extendfrac         = 1.06,
                            dz                 = 2.,
                            ρbg                = 10,
                            nlayers            = 50,
                            ntimesperdecade    = 10,
                            nfreqsperdecade    = 5,
                            Tmax               = -1,
                            nsamples           = -1,
                            nchainsatone       = -1,
                            modelprimary       = false,
                            nchainspersounding = -1) where S<:Sounding

    @assert nsequentialiters  != -1
    @assert nparallelsoundings != -1
    @assert nchainspersounding != -1
    @assert nsamples != - 1
    @assert nchainsatone != -1
    @assert Tmax != -1

    nsoundings = length(soundings)
    opt= deepcopy(opt_in)

    for iter = 1:nsequentialiters
        if iter<nsequentialiters
            ss = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
        else
            ss = (iter-1)*nparallelsoundings+1:nsoundings
        end
        @info "soundings in loop $iter of $nsequentialiters", ss
        r_nothing = Array{Nothing, 1}(undef, length(ss))
        @sync for (i, s) in Iterators.reverse(enumerate(ss))
            pids = (i-1)*nchainspersounding+i:i*(nchainspersounding+1)
            @info "pids in sounding $s:", pids

            aem, = makeoperator(    soundings[s],
                                    zfixed = zfixed,
                                    ρfixed = ρfixed,
                                    zstart = zstart,
                                    extendfrac = extendfrac,
                                    dz = dz,
                                    ρbg = ρbg,
                                    useML = useML,
                                    nlayers = nlayers,
                                    modelprimary = modelprimary,
                                    ntimesperdecade = ntimesperdecade,
                                    nfreqsperdecade = nfreqsperdecade)

            opt = deepcopy(opt_in)
            opt.fdataname = soundings[s].sounding_string*"_"

            @async r_nothing[i] = remotecall_fetch(main, pids[1], opt, aem, collect(pids[2:end]),
                                    Tmax         = Tmax,
                                    nsamples     = nsamples,
                                    nchainsatone = nchainsatone)

        end # @sync
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end

# for gradient based inversion
function loopacrosssoundings(soundings::Array{S, 1}, σstart, σ0; 
                            nsequentialiters   =-1,
                            zfixed             = [-1e5],
                            ρfixed             = [1e12],
                            zstart             = 0.0,
                            extendfrac         = 1.06,
                            dz                 = 2.,
                            ρbg                = 10,
                            nlayers            = 50,
                            ntimesperdecade    = 10,
                            nfreqsperdecade    = 5,
                            modelprimary       = false,
                            regtype            = :R1,
                            nstepsmax          = 10,
                            ntries             = 6,
                            target             = nothing,
                            lo                 = -3.,
                            hi                 = 1.,
                            λ²min              = 0,
                            λ²max              = 8,
                            λ²frac             = 4,
                            β²                 = 0.,
                            ntestdivsλ²        = 50,
                            αmin               = -4, 
                            αmax               = 0, 
                            αfrac              = 4, 
                            ntestdivsα         = 32,
                            regularizeupdate   = false,
                            knownvalue         = 0.7,
                            firstvalue         = :last,
                            κ                  = GP.Mat52(),
                            breakonknown       = true,
                            dobo               = false,
                            ) where S<:Sounding

    @assert nsequentialiters  != -1
    nparallelsoundings = nworkers()
    nsoundings = length(soundings)
    
    for iter = 1:nsequentialiters
        if iter<nsequentialiters
            ss = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
        else
            ss = (iter-1)*nparallelsoundings+1:nsoundings
        end
        @info "soundings in loop $iter of $nsequentialiters", ss
        pids = workers()
        @sync for (i, s) in enumerate(ss)
            aem, = makeoperator(    soundings[s],
                                    zfixed = zfixed,
                                    ρfixed = ρfixed,
                                    zstart = zstart,
                                    extendfrac = extendfrac,
                                    dz = dz,
                                    ρbg = ρbg,
                                    nlayers = nlayers,
                                    modelprimary = modelprimary,
                                    ntimesperdecade = ntimesperdecade,
                                    nfreqsperdecade = nfreqsperdecade,
                                    calcjacobian = true)

            fname = soundings[s].sounding_string*"_gradientinv.dat"
            σstart_, σ0_ = map(x->x*ones(length(aem.ρ)-1), [σstart, σ0])
            @async remotecall_wait(gradientinv, pids[i], σstart_, σ0_, aem,
                                                regtype            = regtype         ,              
                                                nstepsmax          = nstepsmax       ,              
                                                ntries             = ntries          ,              
                                                target             = target          ,              
                                                lo                 = lo              ,              
                                                hi                 = hi              ,              
                                                λ²min              = λ²min           ,              
                                                λ²max              = λ²max           ,              
                                                λ²frac             = λ²frac          ,              
                                                ntestdivsλ²        = ntestdivsλ²     ,              
                                                αmin               = αmin            ,              
                                                αmax               = αmax            ,              
                                                αfrac              = αfrac           ,
                                                β²                 = β²              ,
                                                ntestdivsα         = ntestdivsα      ,              
                                                regularizeupdate   = regularizeupdate,              
                                                knownvalue         = knownvalue      ,              
                                                firstvalue         = firstvalue      ,              
                                                κ                  = κ               ,              
                                                breakonknown       = breakonknown    ,              
                                                dobo               = dobo            ,
                                                fname              = fname           ) 
                

        end # @sync
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end

end
