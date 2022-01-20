module MT1DInversion
using ..AbstractOperator
import ..AbstractOperator.get_misfit
import ..Model, ..Options
using ..MT1D
using Random, PyPlot

mutable struct MT1DZ <: Operator1D
    d_log10_ρ
    d_phase_deg
    σ_log10_ρ
    σ_phase_deg
    freqs
    zboundaries
    irxlayer
end    

function get_misfit(m::Model, opt::Options, F::MT1DZ)
    get_misfit(10 .^m.fstar, opt, F)
end
# above defined function and type signature MUST be defined

function get_misfit(ρ::AbstractArray, opt::Options, F::MT1DZ)
    if !opt.debug
        get_misfit(F.d_log10_ρ, F.d_phase_deg, F.σ_log10_ρ, F.σ_phase_deg, F.freqs, ρ, F.zboundaries, F.irxlayer)
    else
        0.0
    end        
end

function get_misfit(d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg, freqs, ρ, z, irxlayer=1)
    Z = get_Z(freqs, ρ, z, irxlayer)
    ρₐ = MT1D.ρapp(freqs, Z)
    ϕ  = MT1D.phase(Z)
    r = [(d_log10_ρ - log10.(ρₐ))./σ_log10_ρ ; (d_phase_deg - ϕ)./σ_phase_deg]
    r'r/2
end

get_Z(freqs, ρ, z, irxlayer) = MT1D.Z_f(freqs, ρ, diff(z), irxlayer)

function add_noise(ρ, z, freqs; noisefrac=0.05, rseed=1, irxlayer=1)
    Random.seed!(rseed)
    Z = get_Z(freqs, ρ, z, irxlayer)
    σ_real = noisefrac*abs.(Z)/2
    noise = σ_real.*( randn(length(freqs)) + 1im*randn(length(freqs)) )
    Z + noise
end

function create_synthetic(ρ, z, freqs; noisefrac=0.05, rseed=1, showplot=true, logscaledepth=true, showfreq=false, gridalpha=0.5, irxlayer=1)
    Znoisy = add_noise(ρ, z, freqs, noisefrac=noisefrac, rseed=rseed, irxlayer=irxlayer)
    ρₐnoisy = MT1D.ρapp(freqs, Znoisy)
    d_log10_ρ = log10.(ρₐnoisy)   
    d_phase_deg = MT1D.phase(Znoisy)
    σ_log10_ρ = noisefrac/log(10)
    σ_phase_deg = rad2deg(noisefrac/2)
    F = MT1DZ(d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg, freqs, z, irxlayer)
    if showplot
        fig = MT1D.plotmodelcurve(1 ./freqs, ρ, z, logscaledepth=logscaledepth, showfreq=showfreq, irxlayer=irxlayer)
        plotdata(F, fig, iaxis=2, showfreq=showfreq, gridalpha=gridalpha)
    end
    F
end

function plotdata(F::MT1DZ; showfreq=false, gridalpha=0.5, figsize=(6,3))
    fig, ax = plt.subplots(1, 2, sharex=true, figsize=figsize)
    plotdata(F, fig, showfreq=showfreq, gridalpha=gridalpha)
    fig 
end

function plotdata(F::MT1DZ, fig; iaxis=1, showfreq=false, gridalpha=0.5)
    xlabel, abcissa = MT1D.f_or_T(F.freqs, showfreq=showfreq)
    ax = fig.axes
    ax[iaxis].errorbar(abcissa, F.d_log10_ρ, F.σ_log10_ρ, linestyle="none", marker=".", elinewidth=0, capsize=3)
    ax[iaxis+1].errorbar(abcissa, F.d_phase_deg, F.σ_phase_deg, linestyle="none", marker=".", elinewidth=0, capsize=3)
    MT1D.labelaxis(xlabel, ax, iaxis, gridalpha=gridalpha)
    ax[iaxis].set_xscale("log")
    fig.tight_layout()
end    

function plot_posterior(F::MT1DZ, M::AbstractArray; showfreq=false, gridalpha=0.5, logscaledepth=true, 
                            figsize=(10,4), lcolor="nocolor", modelalpha=0.5)
    fig = figure(figsize=(figsize))
    s1 = subplot(131)
    s2 = subplot(132)
    s3 = subplot(133, sharex=s2)
    for m in M
        MT1D.plotmodelcurve(1 ./F.freqs, 10 .^m', F.zboundaries, fig, showfreq=showfreq, irxlayer=F.irxlayer,
                        gridalpha=gridalpha, logscaledepth=logscaledepth, lcolor=lcolor, modelalpha=modelalpha)
    end
    plotdata(F, fig, iaxis=2, showfreq=showfreq, gridalpha=gridalpha)    
end    

end
