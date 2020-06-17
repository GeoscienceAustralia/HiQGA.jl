module CSEM1DInversion
import AbstractOperator.get_misfit
import CSEM1DEr.getfield!
using AbstractOperator, CSEM1DEr
using TransD_GP, PyPlot, LinearAlgebra, CommonToAll

export CSEMRadialEr, makezρ, plotmodelfield!, addnoise

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

function getfield!(m::TransD_GP.Model, csem::CSEMRadialEr)
    copyto!(csem.ρ, csem.nfixed+1:length(csem.ρ), 10 .^m.fstar, 1:length(m.fstar))
    getfield!(csem)
    nothing
end

function get_misfit(m::TransD_GP.Model, opt::TransD_GP.Options, csem::CSEMRadialEr)
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
                        ;figsize=(10,5))
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
                  )
    @assert all((nfixed, dz, extendfrac) .> 0)
    CSEM1DEr.getfield!(F, z, ρ)
    d = F.Er + sqrt(2)*noisefrac*abs.(F.Er).*randn(ComplexF64, size(F.Er))
    d[abs.(d).<noisefloor] .= NaN
    plotmodelfield!(F, z, ρ, d, figsize=figsize, nfixed=nfixed,
                                dz=dz, extendfrac=extendfrac)
    return d, sqrt(2)*noisefrac*abs.(d)
end

function plotmodelfield!(F::CSEM1DEr.RadialEr, z, ρ, d::Array{ComplexF64}, σ::Array{Float64};
                        figsize=(12,4), nfixed=-1, dz=-1., extendfrac=-1.)
    @assert all((nfixed, dz, extendfrac) .> 0)
    f, ax = plt.subplots(1, 3, figsize=figsize)
    ax[1].step(log10.(ρ[2:end]), z[2:end])
    if dz > 0
        axn = ax[1].twinx()
        ax[1].get_shared_y_axes().join(ax[1],axn)
        axn.step(log10.(ρ[2:end]), z[2:end])
        yt = ax[1].get_yticks()[ax[1].get_yticks().>0]
        axn.set_yticks(yt)
        axn.set_ylim(ax[1].get_ylim()[end:-1:1])
        axn.set_yticklabels(string.(Int.(round.(getn.(yt .- z[nfixed+1], dz, extendfrac)))))
    end
    CSEM1DEr.getfield!(F, z, ρ)
    ax[2].semilogy(F.rRx, abs.(F.Er))
    ax[3].plot(F.rRx, unwrap(angle.(F.Er)))
    for lfreq in eachindex(F.freqs)
        ax[2].errorbar(F.rRx, abs.(d[:,lfreq]), yerr = 2*0.707abs.(σ[:,lfreq]),
                        linestyle="none", marker=".", elinewidth=0, capsize=3)
        ax[3].errorbar(F.rRx, unwrap(angle.(d[:,lfreq])),
                                yerr = 2*0.707abs.(σ[:,lfreq])./abs.(d[:,lfreq]),
                        linestyle="none", marker=".", elinewidth=0, capsize=3)
    end
    ax[1].grid()
    ax[2].grid()
    ax[3].grid()
    nicenup(f)
end

end
