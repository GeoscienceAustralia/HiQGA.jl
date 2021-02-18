module SkyTEM1DInversion
import ..AbstractOperator.get_misfit
# include("AEM_VMD_HMD.jl")
using ..AbstractOperator, ..AEM_VMD_HMD, Statistics
using PyPlot, LinearAlgebra, ..CommonToAll, MAT, Random, DelimitedFiles

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
end

function dBzdt(Flow       :: AEM_VMD_HMD.HField,
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
    dBzdt(dlow, dhigh, useML, σlow, σhigh,
    Flow, Fhigh, z, nfixed, copy(ρ), selectlow, selecthigh, ndatalow, ndatahigh)
end

mutable struct SkyTEMsoundingData
    sounding_string :: String
    X :: Float64
    Y :: Float64
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
                            sounding_string="sounding", X=nothing, Y=nothing,
                            linenum=nothing, fid=nothing)
    @assert rRx > 0 && rTx > 0
    @assert zRxLM <0 && zTxLM <0
    @assert zRxHM <0 && zTxHM <0
    @assert all(lowpassfcs .> 0)
    @assert all(diff(LM_times) .>0 )
    @assert all(diff(HM_times) .>0 )
    @assert all(diff(LM_ramp[:,1]) .>0 )
    @assert all(diff(HM_ramp[:,1]) .>0 )
    @assert all(LM_noise .>0 )
    @assert all(HM_noise .>0 )
    @assert length(LM_data) == length(LM_noise)
    @assert length(HM_data) == length(HM_noise)
    SkyTEMsoundingData(sounding_string, X, Y, fid, linenum, rRx, zRxLM, zTxLM, zRxHM, zTxHM, rTx,
    lowpassfcs, LM_times, LM_ramp, HM_times, HM_ramp, LM_noise, HM_noise,
    LM_data, HM_data)
    # @show (sounding_string, rRx, zRxLM, zTxLM, zRxHM, zTxHM, rTx,
    # lowpassfcs, LM_times, LM_ramp, HM_times, HM_ramp, LM_noise, HM_noise,
    # LM_data, HM_data)
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
    units=1e-12,
    figsize = (9,7),
    makesounding = false,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.03,
    X = -1,
    Y = -1,
    fid = -1,
    linenum = -1)
    @assert frame_height > 0
    @assert frame_dz > 0
    @assert frame_dx > 0
    @assert frame_dy > 0
    @assert all(LM_Z .> 0)
    @assert all(HM_Z .> 0)
    @assert X > 0
    @assert Y > 0
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
    fiducial = soundings[:,fid]
    whichline = soundings[:,linenum]
    d_LM = soundings[:,LM_Z[1]:LM_Z[2]]
    d_HM = soundings[:,HM_Z[1]:HM_Z[2]]
    zTx = soundings[:,frame_height]
    zRx = -(zTx + soundings[:,frame_dz])
    zTx = -zTx
    rRx = sqrt.(soundings[:,frame_dx].^2 + soundings[:,frame_dy].^2)

    @info "reading $fname_specs_halt"
    include(fname_specs_halt)
    @assert size(d_LM, 2) == length(LM_times)
    @assert size(d_HM, 2) == length(HM_times)
    @assert size(d_LM, 2) == length(LM_noise)
    @assert size(d_HM, 2) == length(HM_noise)
    LM_noise[:] .*= units
    HM_noise[:] .*= units
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
            dₗₘ, dₕₘ = vec(d_LM[is,:]), vec(d_HM[is,:])
            σLM = sqrt.((multnoise*dₗₘ).^2 + LM_noise.^2)
            σHM = sqrt.((multnoise*dₕₘ).^2 + HM_noise.^2)
            s_array[is] = SkyTEMsoundingData(rRx=rRx[is], zRxLM=zRx[is], zTxLM=zTx[is],
                zRxHM=zRx[is], zTxHM=zTx[is], rTx=rTx, lowpassfcs=lowpassfcs,
                LM_times=LM_times, LM_ramp=LM_ramp,
                HM_times=HM_times, HM_ramp=HM_ramp,
                LM_noise=σLM, HM_noise=σHM, LM_data=dₗₘ, HM_data=dₕₘ,
                sounding_string="sounding_$(l)_$f",
                X=easting[is], Y=northing[is], fid=f,
                linenum=l)
        end
        return s_array
    end
end

function getfield!(m::Model, aem::dBzdt)
    copyto!(aem.ρ, aem.nfixed+1:length(aem.ρ), 10 .^m.fstar, 1:length(m.fstar))
    AEM_VMD_HMD.getfieldTD!(aem.Flow,  aem.z, aem.ρ)
    AEM_VMD_HMD.getfieldTD!(aem.Fhigh, aem.z, aem.ρ)
    nothing
end

function getfield!(m::Array{Float64}, aem::dBzdt)
    copyto!(aem.ρ, aem.nfixed+1:length(aem.ρ), 10 .^m, 1:length(m))
    AEM_VMD_HMD.getfieldTD!(aem.Flow,  aem.z, aem.ρ)
    AEM_VMD_HMD.getfieldTD!(aem.Fhigh, aem.z, aem.ρ)
    nothing
end

function get_misfit(m::Model, opt::Options, aem::dBzdt)
    chi2by2 = 0.0
    if !opt.debug
        getfield!(m, aem)
        chi2by2 += getchi2by2(aem.Flow.dBzdt, aem.dlow,
                    aem.σlow, aem.selectlow, aem.useML, aem.ndatalow)
        chi2by2 += getchi2by2(aem.Fhigh.dBzdt, aem.dhigh,
                    aem.σhigh, aem.selecthigh, aem.useML, aem.ndatahigh)
    end
    return chi2by2
end

function getchi2by2(dBzdt, d, σ, select, useML, ndata)
    r, d, s, idx = dBzdt, d, σ, select
    r .= (r - d)./s
    if useML
        chi2by2 = 0.5*ndata[ifreq]*log(norm(r[idx])^2)
    else
        chi2by2 = 0.5*norm(r[idx])^2
    end
end

function plotmodelfield_skytem!(Flow::AEM_VMD_HMD.HField, Fhigh::AEM_VMD_HMD.HField,
                         z::Array{Float64, 1}, ρ::Array{Float64, 1}
                        ;figsize=(10,5))
    f, ax = plt.subplots(1, 2, figsize=figsize)
    ax[1].step(log10.(ρ[2:end]), z[2:end])
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
                  noisefrac  = 0.05,
                  noisefloorlow = μ₀*1e-14,
                  noisefloorhigh = μ₀*1e-14,
                  dz = -1.,
                  extendfrac = -1.,
                  nfixed = -1,
                  figsize=(12,4)
                  )
    @assert all((nfixed, dz, extendfrac) .> 0)
    AEM_VMD_HMD.getfieldTD!(Flow, z, ρ)
    AEM_VMD_HMD.getfieldTD!(Fhigh, z, ρ)
    dlow  = Flow.dBzdt + noisefrac*abs.(Flow.dBzdt).*randn(size(Flow.dBzdt))
    dhigh = Fhigh.dBzdt + noisefrac*abs.(Fhigh.dBzdt).*randn(size(Fhigh.dBzdt))
    dlow[abs.(dlow).<noisefloorlow] .= NaN
    dhigh[abs.(dhigh).<noisefloorhigh] .= NaN
    σlow = noisefrac*abs.(dlow)
    σhigh = noisefrac*abs.(dhigh)
    plotmodelfield!(Flow, Fhigh, z, ρ, dlow, dhigh, σlow, σhigh,
                                figsize=figsize, nfixed=nfixed,
                                dz=dz, extendfrac=extendfrac)
    # returned data is dBzdt not H if there is a μ multiplied
    return μ₀.*(dlow, dhigh, σlow, σhigh)
end

function plotmodelfield!(Flow::AEM_VMD_HMD.HField, Fhigh::AEM_VMD_HMD.HField,
                        z, ρ, dlow, dhigh, σlow, σhigh;
                        figsize=(12,4), nfixed=-1, dz=-1., extendfrac=-1.)
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
    AEM_VMD_HMD.getfieldTD!(Flow, z, ρ)
    AEM_VMD_HMD.getfieldTD!(Fhigh, z, ρ)
    ax[2].loglog(Flow.times,μ₀*Flow.dBzdt, label="low moment")
    ax[2].loglog(Fhigh.times,μ₀*Fhigh.dBzdt, label="high moment")
    ax[2].errorbar(Flow.times, μ₀*vec(dlow), yerr = μ₀*2abs.(vec(σlow)),
                        linestyle="none", marker=".", elinewidth=0, capsize=3)
    ax[2].errorbar(Fhigh.times, μ₀*vec(dhigh), yerr = μ₀*2abs.(vec(σhigh)),
                        linestyle="none", marker=".", elinewidth=0, capsize=3)
    ax[1].grid()
    ax[2].grid()
    nicenup(f)
end

function plotmodelfield!(aem::dBzdt, Ρ::Vector{Array{Float64}};
                        figsize=(8,5), dz=-1., onesigma=true,
                        extendfrac=-1., fsize=10, alpha=0.1)
    @assert all((dz, extendfrac) .> 0)
    sigma = onesigma ? 1.0 : 2.0
    f = figure(figsize=figsize)
    ax = Vector{PyPlot.PyObject}(undef, 3)
    ax[1] = subplot(121)
    ρmin, ρmax = extrema(vcat(Ρ...))
    delρ = ρmax - ρmin
    ax[1].set_xlim(ρmin-0.1delρ,ρmax+0.1delρ)
    nfixed, z = aem.nfixed, aem.z
    ax[1].plot([ρmin-0.1delρ,ρmax+0.1delρ], z[nfixed+1]*[1., 1], color="b")
    ax[2] = subplot(122)
    Flow = aem.Flow
    dlow, σlow = aem.dlow, aem.σlow
    Fhigh = aem.Fhigh
    dhigh, σhigh = aem.dhigh, aem.σhigh
    ax[2].errorbar(Flow.times, μ₀*dlow, yerr = μ₀*sigma*abs.(σlow),
                        linestyle="none", marker=".", elinewidth=0, capsize=3, label="low moment")
    ax[2].errorbar(Fhigh.times, μ₀*dhigh, yerr = μ₀*sigma*abs.(σhigh),
                        linestyle="none", marker=".", elinewidth=0, capsize=3, label="high moment")
    for ρ in Ρ
        getfield!(ρ,  aem)
        Flow.dBzdt[.!aem.selectlow] .= NaN
        Fhigh.dBzdt[.!aem.selecthigh] .= NaN
        ax[1].step(log10.(aem.ρ[2:end]), aem.z[2:end], "-k", alpha=alpha)
        ax[2].loglog(Flow.times,μ₀*Flow.dBzdt, "k", alpha=alpha, markersize=2)
        ax[2].loglog(Fhigh.times,μ₀*Fhigh.dBzdt, "k", alpha=alpha, markersize=2)
    end
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
    axn.set_ylabel("Depth index", rotation=-90)
    ax[1].set_xlabel("Log₁₀ρ")
    ax[1].set_title("Model")
    ax[2].set_ylabel(L"dBzdt \; V/(A.m^4)")
    ax[2].set_xlabel("Time (s)")
    ax[2].set_title("Transient response")
    ax[2].legend()
    ax[2].grid()
    ax[1].invert_xaxis()
    nicenup(f, fsize=fsize)
end

function makeoperator(fdataname::String;
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
                       plotfield = false)
    @assert extendfrac > 1.0
    @assert dz > 0.0
    @assert ρbg > 0.0
    @assert nlayers > 1
    nmax = nlayers+1

    zall, znall, zboundaries = setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=showgeomplot)
    z, ρ, nfixed = makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
    ρ[z.>=zstart] .= ρbg
    ##  geometry and modeling parameters
    file = matopen(fdataname)
    rRx = read(file, "rRx")
    zRxLM = read(file, "LM_zRx")
    zTxLM = read(file, "LM_zTx")
    zRxHM = read(file, "HM_zRx")
    zTxHM = read(file, "HM_zTx")
    rTx = read(file, "rTxLoop")
    lowpassfcs = read(file, "lowPassFilters")
    # Note that the receiver depth needs to be in same model layer as transmitter.
    ## LM times and ramp
    LM_times = read(file, "LM_times")[:]
    LM_ramp = read(file, "LM_ramp")
    ## HM times and ramp
    HM_times = read(file, "HM_times")[:]
    HM_ramp = read(file, "HM_ramp")
    ## LM operator
    Flm = AEM_VMD_HMD.HFieldDHT(
                          ntimesperdecade = ntimesperdecade,
                          nfreqsperdecade = nfreqsperdecade,
                          lowpassfcs = lowpassfcs,
                          times  = LM_times,
                          ramp   = LM_ramp,
                          nmax   = nmax,
                          zTx    = zTxLM,
                          rRx    = rRx,
                          rTx    = rTx,
                          zRx    = zRxLM)
    ## HM operator
    Fhm = AEM_VMD_HMD.HFieldDHT(
                          ntimesperdecade = ntimesperdecade,
                          nfreqsperdecade = nfreqsperdecade,
                          lowpassfcs = lowpassfcs,
                          times  = HM_times,
                          ramp   = HM_ramp,
                          nmax   = nmax,
                          zTx    = zTxHM,
                          rRx    = rRx,
                          rTx    = rTx,
                          zRx    = zRxHM)
    ## data and high altitude noise
    LM_data = read(file, "d_LM")
    HM_data = read(file, "d_HM")
    LM_noise = read(file, "sd_LM")
    HM_noise = read(file, "sd_HM")
    ## create operator
    dlow, dhigh, σlow, σhigh = (LM_data, HM_data, LM_noise, HM_noise)./μ₀
    aem = dBzdt(Flm, Fhm, vec(dlow), vec(dhigh),
                                      vec(σlow), vec(σhigh), z=z, ρ=ρ, nfixed=nfixed)
    plotfield && plotmodelfield!(Flm, Fhm, z, ρ, dlow, dhigh, σlow, σhigh;
                          figsize=(12,4), nfixed=nfixed, dz=dz, extendfrac=extendfrac)
    aem, znall
end

function makeoperator( sounding::SkyTEMsoundingData;
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
                          zRx    = sounding.zRxLM)
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
                          zRx    = sounding.zRxHM)
    ## create operator
    dlow, dhigh, σlow, σhigh = (sounding.LM_data, sounding.HM_data, sounding.LM_noise, sounding.HM_noise)./μ₀
    aem = dBzdt(Flm, Fhm, vec(dlow), vec(dhigh),
                                      vec(σlow), vec(σhigh), z=z, ρ=ρ, nfixed=nfixed)
    plotfield && plotmodelfield!(Flm, Fhm, z, ρ, dlow, dhigh, σlow, σhigh;
                          figsize=(12,4), nfixed=nfixed, dz=dz, extendfrac=extendfrac)
    aem, znall
end

function make_tdgp_statmode_opt(;
                    rseed = nothing,
                    znall = znall,
                    fileprefix = "sounding",
                    nmin = 2,
                    nmax = 40,
                    K = GP.Mat32(),
                    demean = true,
                    sdpos = 0.05,
                    sdprop = 0.05,
                    fbounds = [-0.5 2.5],
                    λ = [2],
                    δ = 0.1,
                    pnorm = 2,
                    save_freq = 25
                    )
    sdev_pos = [sdpos*abs(diff([extrema(znall)...])[1])]
    sdev_prop = sdprop*diff(fbounds, dims=2)[:]
    xall = permutedims(collect(znall))
    xbounds = permutedims([extrema(znall)...])

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
                            pnorm = pnorm,
                            quasimultid = false,
                            K = K,
                            save_freq = save_freq,
                            needλ²fromlog = needλ²fromlog,
                            updatenonstat = updatenonstat,
                            dispstatstoscreen = false
                            )

    ## Initialize options for the dummy nonstationary properties GP
    optdummy = OptionsNonstat(opt,
                            nmin = 2,
                            nmax = 3,
                            fbounds = fbounds,
                            δ = δ,
                            demean = demean,
                            sdev_prop = sdev_prop,
                            sdev_pos = sdev_pos,
                            pnorm = pnorm,
                            K = K
                            )

    opt, optdummy

end


end