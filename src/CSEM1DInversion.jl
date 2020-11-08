module CSEM1DInversion
import ..AbstractOperator.get_misfit
# include("CSEM1DEr.jl")
import ..CSEM1DEr.getfield!
using ..AbstractOperator, ..CSEM1DEr
using PyPlot, LinearAlgebra, ..CommonToAll

import ..Model, ..Options

export CSEMRadialEr, plotmodelfield!, addnoise

mutable struct CSEMRadialEr<:Operator1D
    d      :: Array{ComplexF64}
    useML  :: Bool
    σ      :: Array{Float64}
    F      :: CSEM1DEr.RadialEr
    z      :: Array{Float64, 1}
    nfixed :: Int
    ρ      :: Array{Float64, 1}
    select :: Array{Bool}
    ndata  :: Array{Int, 1}
end

function makezρ(zboundaries::Array{Float64, 1};
          zfixed      = [-1e6,    0.],
          ρfixed      = [1e12     0.3])
    @assert length(zfixed) == length(ρfixed)
    z = [zfixed..., zboundaries...]
    @assert all(diff(z).>0)
    ρ = zeros(length(z))
    nfixed = length(ρfixed)
    ρ[1:nfixed] .= ρfixed
    z, ρ, nfixed
end

function CSEMRadialEr(F          :: CSEM1DEr.RadialEr,
                      d          :: Array{ComplexF64},
                      σ          :: Array{Float64};
                      useML  = false,
                      z = [-1.],
                      ρ = [-1.],
                      nfixed = -1
                    )
    @assert length(F.thickness) >= length(z)
    @assert size(σ) == size(d)
    ndata = vec(sum(.!isnan.(d), dims=1))
    select = .!isnan.(d)
    CSEMRadialEr(d, useML, σ, F, z, nfixed, copy(ρ), select, ndata)
end

getfield!(csem::CSEMRadialEr) = getfield!(csem.F, csem.z, csem.ρ)

function getfield!(m::Model, csem::CSEMRadialEr)
    copyto!(csem.ρ, csem.nfixed+1:length(csem.ρ), 10 .^m.fstar, 1:length(m.fstar))
    getfield!(csem)
    nothing
end

function getfield!(m::Array{Float64}, csem::CSEMRadialEr)
    copyto!(csem.ρ, csem.nfixed+1:length(csem.ρ), 10 .^m, 1:length(m))
    getfield!(csem)
    nothing
end

function get_misfit(m::Model, opt::Options, csem::CSEMRadialEr)
    chi2by2 = 0.0
    if !opt.debug
        ndata = csem.ndata
        getfield!(m, csem)
            @views begin
            for ifreq in eachindex(csem.F.freqs)
                r, d, s, idx = csem.F.Er[:,ifreq], csem.d[:,ifreq], csem.σ[:,ifreq], csem.select[:,ifreq]
                r .= (r - d)./s
                if csem.useML
                    chi2by2 += ndata[ifreq]*log(norm(r[idx])^2)
                else
                    chi2by2 += norm(r[idx])^2
                end
            end
        end
    end
    return chi2by2
end

function plotmodelfield!(F::CSEM1DEr.RadialEr, z::Array{Float64, 1}, ρ::Array{Float64, 1}
                        ;figsize=(12,4))
    f, ax = plt.subplots(1, 3, figsize=figsize)
    ax[1].step(log10.(ρ[2:end]), z[2:end])
    CSEM1DEr.getfield!(F, z, ρ)
    ax[2].semilogy(F.rRx, abs.(F.Er))
    ax[3].plot(F.rRx, unwrap(angle.(F.Er)))
    ax[1].grid()
    ax[1].invert_yaxis()
    ax[2].grid()
    ax[3].grid()
    nicenup(f)
end

function addnoise(F::CSEM1DEr.RadialEr, z::Array{Float64, 1}, ρ::Array{Float64, 1};
                  noisefrac  = 0.05,
                  noisefloor = 1e-14,
                  dz = -1.,
                  extendfrac = -1.,
                  nfixed = -1,
                  figsize=(8,3)
                  )
    @assert all((nfixed, dz, extendfrac) .> 0)
    CSEM1DEr.getfield!(F, z, ρ)
    d = F.Er + sqrt(2)*noisefrac*abs.(F.Er).*randn(ComplexF64, size(F.Er))
    d[abs.(d).<noisefloor] .= NaN
    σ = sqrt(2)*noisefrac*abs.(d)
    plotmodelfield!(F, z, ρ, d, σ, figsize=figsize, nfixed=nfixed,
                                dz=dz, extendfrac=extendfrac)
    return d, σ
end

function plotmodelfield!(F::CSEM1DEr.RadialEr, z, ρ::Vector{Float64}, d::Array{ComplexF64}, σ::Array{Float64};
                        figsize=(8,3), nfixed=-1, dz=-1., extendfrac=-1., fsize=8, onesigma=true)
    @assert all((nfixed, dz, extendfrac) .> 0)
    sigma = onesigma ? 1.0 : 2.0
    f = figure(figsize=figsize)
    ax = Vector{PyPlot.PyObject}(undef, 3)
    ax[1] = subplot(131)
    ax[1].step(log10.(ρ[2:end]), z[2:end])
    ρmin, ρmax = extrema(log10.(ρ[2:end]))
    delρ = ρmax - ρmin
    ax[1].set_xlim(ρmin-0.1delρ,ρmax+0.1delρ)
    ax[1].plot([ρmin-0.1delρ,ρmax+0.1delρ], z[nfixed+1]*[1., 1])
    if dz > 0
        axn = ax[1].twinx()
        ax[1].get_shared_y_axes().join(ax[1],axn)
        yt = ax[1].get_yticks()[ax[1].get_yticks().>=z[nfixed+1]]
        axn.set_yticks(yt)
        axn.set_ylim(ax[1].get_ylim()[end:-1:1])
        axn.set_yticklabels(string.(Int.(round.(getn.(yt .- z[nfixed+1], dz, extendfrac)))))
    end
    ax[2] = subplot(132)
    ax[3] = subplot(133, sharex=ax[2], sharey=ax[2])
    for lfreq in eachindex(F.freqs)
        yerr = sigma*0.707abs.(σ[:,lfreq])
        ax[2].errorbar(F.rRx, abs.(real(d[:,lfreq])), yerr = yerr,
                        linestyle="none", marker=".", elinewidth=0, capsize=3, alpha=0.5, markersize=5)
        ax[3].errorbar(F.rRx, abs.(imag(d[:,lfreq])), yerr = yerr,
                        linestyle="none", marker=".", elinewidth=0, capsize=3, alpha=0.5, markersize=5)
    end
    CSEM1DEr.getfield!(F, z, ρ)
    ax[2].semilogy(F.rRx, abs.(real(F.Er)), "-k", alpha=0.5)
    ax[3].semilogy(F.rRx, abs.(imag(F.Er)), "-k", alpha=0.5)
    ax[1].grid()
    ax[1].set_ylabel("Depth m")
    axn.set_ylabel("Depth index", rotation=-90)
    ax[1].set_xlabel("Log₁₀ρ")
    ax[1].set_title("Model")
    ax[2].set_ylabel(L"E_r \; V/(A.m^2)")
    ax[2].set_xlabel("Receiver range m")
    ax[2].set_title("Real "*L"E_r")
    ax[2].grid()
    ax[3].grid()
    ax[3].set_title("Imaginary "*L"E_r")
    ax[3].set_ylabel(L"E_r \; V/(A.m^2)")
    ax[3].set_xlabel("Receiver range m")
    ampmin = minimum(abs.([real(d[.!isnan.(d)])...,imag(d[.!isnan.(d)])...]))
    yl = ax[2].get_ylim()
    yl = (0.9*ampmin, yl[2])
    ax[2].set_ylim(yl)
    nicenup(f, fsize=fsize)
    box = ax[3].get_position()
    box.x0 = box.x0 - 0.05
    box.x1 = box.x1 - 0.05
    ax[3].set_position(box)
end

function plotmodelfield!(csem::CSEMRadialEr, Ρ::Vector{Array{Float64}};
                        figsize=(8,3), dz=-1., onesigma=true,
                        extendfrac=-1., fsize=8, alpha=0.1)
    @assert all((dz, extendfrac) .> 0)
    sigma = onesigma ? 1.0 : 2.0
    f = figure(figsize=figsize)
    ax = Vector{PyPlot.PyObject}(undef, 3)
    ax[1] = subplot(131)
    ρmin, ρmax = extrema(vcat(Ρ...))
    delρ = ρmax - ρmin
    ax[1].set_xlim(ρmin-0.1delρ,ρmax+0.1delρ)
    nfixed, z = csem.nfixed, csem.z
    ax[1].plot([ρmin-0.1delρ,ρmax+0.1delρ], z[nfixed+1]*[1., 1], color="b")
    ax[2] = subplot(132)
    ax[3] = subplot(133, sharex=ax[2], sharey=ax[2])
    F = csem.F
    d, σ = csem.d, csem.σ
    for lfreq in eachindex(F.freqs)
        yerr = sigma*0.707abs.(σ[:,lfreq])
        ax[2].errorbar(F.rRx, abs.(real(d[:,lfreq])), yerr = yerr, label=string(F.freqs[lfreq])*" Hz",
                        linestyle="none", marker=".", elinewidth=0, capsize=3, alpha=0.3, markersize=5)
        ax[3].errorbar(F.rRx, abs.(imag(d[:,lfreq])), yerr = yerr,
                        linestyle="none", marker=".", elinewidth=0, capsize=3, alpha=0.3, markersize=5)
    end
    for ρ in Ρ
        getfield!(ρ, csem)
        F.Er[.!csem.select] .= NaN
        ax[1].step(log10.(csem.ρ[2:end]), csem.z[2:end], "-k", alpha=alpha)
        ax[2].semilogy(F.rRx, abs.(real(F.Er)), ".k", alpha=alpha, markersize=2)
        ax[3].semilogy(F.rRx, abs.(imag(F.Er)), ".k", alpha=alpha, markersize=2)
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
    ax[2].set_ylabel(L"E_r \; V/(A.m^2)")
    ax[2].set_xlabel("Receiver range m")
    ax[2].set_title("Real "*L"E_r")
    ax[2].legend()
    ax[2].grid()
    ax[3].grid()
    ax[3].set_title("Imaginary "*L"E_r")
    ax[3].set_ylabel(L"E_r \; V/(A.m^2)")
    ax[3].set_xlabel("Receiver range m")
    ampmin = minimum(abs.([real(d[.!isnan.(d)])...,imag(d[.!isnan.(d)])...]))
    yl = ax[2].get_ylim()
    yl = (0.9*ampmin, yl[2])
    ax[2].set_ylim(yl)
    nicenup(f, fsize=fsize)
    box = ax[3].get_position()
    box.x0 = box.x0 - 0.05
    box.x1 = box.x1 - 0.05
    ax[3].set_position(box)
end

end
