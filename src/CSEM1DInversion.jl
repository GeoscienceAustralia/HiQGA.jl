module CSEM1DInversion
import AbstractOperator.get_misfit
using AbstractOperator, CSEM1DEr
using TransD_GP, PyPlot, LinearAlgebra, CommonToAll

export CSEMRadialEr

mutable struct CSEMRadialEr<:Operator1D
    d      :: Array{Float64}
    useML  :: Bool
    σ      :: Float64
    F      :: RadialEr
    r      :: Array{ComplexF64}
    d      :: Array{ComplexF64}
    z      :: Array{Float64, 1}
    ρ      :: Array{Float64, 1}
end

function CSEMRadialEr(d::Array{Float64}, zfixed, zsetup, ρfixed;
                        zTx = 975.
                        zRx = 1000.,
                        rRx = collect(500:32:5000)
                        RxAzim = 0.,
                        TxDip = 0.,
                        nmax = 200,
                        useML=false,
                        σ=1.0)
    F = CSEM1DEr.RadialErLagged(  zTx    = [zTx],
                                  rRx    = rRx,
                                  freqs  = freqs,
                                  zRx    = [zRx],
                                  RxAzim = RxAzim,
                                  TxDip  = TxDip,
                                  nmax   = nmax)
    z = [zfixed, zsetup]
    ρ = [ρfixed, zeros(nmax)...]
    CSEMRadialEr(d, useML, σ, F, zeros(ComplexF64, size(F.Er), z, ρ)
end

function get_misfit(m::TransD_GP.Model, opt::TransD_GP.Options, csem::CSEMRadialEr)
    chi2by2 = 0.0
    if !opt.debug
        d, Er = csem.d, csem.F.Er
        select = .!isnan.(d)
        getfield!(csem.F, z::Array{Float64, 1}, ρ::Array{Float64, 1})
        csem.r = Er[select] - d[select]
        if line.useML
            N = sum(select)
            chi2by2 = 0.5*N*log(norm(r)^2)
        else
            chi2by2 = r'*r/(2line.σ^2)
        end
    end
    return chi2by2
end

end
