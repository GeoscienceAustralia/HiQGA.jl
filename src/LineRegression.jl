module LineRegression
import AbstractOperator.get_misfit
using AbstractOperator
using TransD_GP, PyPlot, LinearAlgebra

export Line, linetestfunction

mutable struct Line<:Operator1D
    d     :: Array{Float64}
    useML :: Bool
    σ     :: Float64
end

function Line(d::Array{Float64, 1} ;useML=false, σ=1.0)
    Line(d, useML, σ)
end

function get_misfit(m::TransD_GP.Model, opt::TransD_GP.Options, line::Line)
    chi2by2 = 0.0
    if !opt.debug
        d = line.d
        select = .!isnan.(d[:])
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

function linetestfunction(;c=0.25,ngrid=200)
    xx = LinRange(-1,1,ngrid)
    y = zeros(size(xx))
    for (i, x) in enumerate(xx)
        y[i] = x <= -c ? -1 - 2(x+c)*(x+c) : 2 + 2*x*x
    end
    xx, y
end

end
