module CSEM1DInversion
import AbstractOperator.get_misfit
using AbstractOperator
import CSEM1DEr.RadialEr
using TransD_GP, PyPlot, LinearAlgebra, CommonToAll

export CSEMRadialEr

mutable struct CSEMRadialEr<:Operator1D
    d     :: Array{Float64}
    useML :: Bool
    σ     :: Float64
    F     :: RadialEr
    r     :: Array{ComplexF64}
end

function CSEMRadialEr(d::Array{Float64}, F::RadialEr; useML=false, σ=1.0,)
    CSEMRadialEr(d, useML, σ, F, zeros{ComplexF64, size(F.d)})
end

function geomprogdepth(n, dy, c)
    dy*(1.0-c^n)/(1-c)
end

function getn(z, dy, c)
    log(1 - z/dy*(1-c))/log(c)
end

function setupz(;n=100, dz=10, extendfrac=1.01)
    znrange            = 1.0:n
    zboundaries        = geomprogdepth.(znrange, dz, extendfrac)
    thickness          = [zboundaries[1]; diff(zboundaries)[1:end-1]]
    zall               = [zboundaries[1]/2; 0.5*(zboundaries[1:end-1] + zboundaries[2:end])]
    znall              = getn.(zall, dz, extendfrac)

    figure()
    plot(znall, zall)
    xlabel("depth index")
    ylabel("depth associated")
    grid()
    nicenup(gcf())

    f, ax = plt.subplots(1, 2, figsize=(10,5))
    ax[1].stem(zboundaries[1:end-1], zboundaries[1:end-1], markerfmt="")
    ax[1].stem(zall, zall, "k--", markerfmt=" ")
    ax[1].set_xlabel("depth m")
    ax[1].set_ylabel("depth m")

    ax[2].stem(znrange[1:end-1], znrange[1:end-1], markerfmt="")
    ax[2].stem(znall, znall, "k--", markerfmt=" ")
    ax[2].set_ylabel("depth index")
    ax[2].set_xlabel("depth index")
    nicenup(gcf())

    f, ax = plt.subplots(1, 2, figsize=(10,5), sharey=true)
    ax[1].stem(zall[1:end-1],thickness, "k--", markerfmt=" ")
    ax[1].set_ylabel("thickness m")
    ax[1].yaxis.grid(which="major")
    ax[2].stem(znall[1:end-1],thickness, "k--", markerfmt=" ")
    ax[2].set_xlabel("depth index")
    ax[2].yaxis.grid(which="major")
    nicenup(f)
end

function get_misfit(m::TransD_GP.Model, opt::TransD_GP.Options, csem::CSEMRadialEr)
    chi2by2 = 0.0
    if !opt.debug
        d = csem.d
        select = .!isnan.(d[:])
        getfield!(csem.F, z::Array{Float64, 1}, ρ::Array{Float64, 1})
        r = m.fstar[:][select] - d[select]
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
