module CSEM1DInversion
import AbstractOperator.get_misfit
using AbstractOperator, CSEM1DEr
using TransD_GP, PyPlot, LinearAlgebra, CommonToAll

export CSEMRadialEr

mutable struct CSEMRadialEr<:Operator1D
    d      :: Array{ComplexF64}
    useML  :: Bool
    σ      :: Array{Float64}
    F      :: RadialEr
    z      :: Array{Float64, 1}
    nfixed :: Int
    ρ      :: Array{Float64, 1}
    select :: Array{Bool}
    ndata  :: Array{Int, 1}
end

function CSEMRadialEr(d::Array{Float64},
                    zfixed::Array{Float64, 1},
                    ρfixed::Array{Float64, 1};
                    zTx = 975.,
                    zRx = 1000.,
                    rRx = collect(500:32:5000),
                    RxAzim = 0.,
                    TxDip = 0.,
                    nmax = 200,
                    useML=false,
                    σ=[1.0])
    @assert length(zfixed) == length(ρfixed)
    z = [zfixed, zsetup...]
    @assert all(diff(z).>0)
    @assert nmax >= length(z)
    @assert size(σ) == size(d)
    ρ = zeros(length(z))
    nfixed = length(ρfixed)
    ndata = permutedims(sum(.!isnan.(d[:,ifreq])))
    select = .!isnan(d)
    ρ[1:nfixed] .= ρfixed
    F = CSEM1DEr.RadialErLagged(  zTx    = [zTx],
                                  rRx    = rRx,
                                  freqs  = freqs,
                                  zRx    = [zRx],
                                  RxAzim = RxAzim,
                                  TxDip  = TxDip,
                                  nmax   = nmax)
    CSEMRadialEr(d, useML, σ, F, z, nfixed, ρ, select, ndata)
end

function get_misfit(m::TransD_GP.Model, opt::TransD_GP.Options, csem::CSEMRadialEr)
    chi2by2 = 0.0
    if !opt.debug
        csem.ρ[csem.nfixed+1:end] .= 10 .^m.fstar
        d, Er, select, ndata, σ = csem.d, csem.F.Er, csem.select, csem.ndata, csem.σ
        getfield!(csem.F, csem.z, csem.ρ)
        @views for ifreq in eachindex(csem.F.freqs)
            r, d, s, idx = Er[:,ifreq], d[:,ifreq], σ[:,ifreq], select[:,ifreq]
            r .= (r - d)./s
            if csem.useML
                chi2by2 += ndata[ifreq]*log(norm(r[idx])^2)
            else
                chi2by2 += norm(r[idx])^2
            end
        end
    end
    return chi2by2
end

end
