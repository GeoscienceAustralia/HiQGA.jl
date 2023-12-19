module LineRegression
import ..AbstractOperator.get_misfit
using ..AbstractOperator, ..CommonToAll
using PyPlot, LinearAlgebra, StatsBase, SparseArrays

import ..Model, ..Options
import ..AbstractOperator.getresidual

export Line, makehist

mutable struct Line<:Operator1D
    d      :: Array{Float64, 1}
    useML  :: Bool
    σ      :: Union{Array{Float64, 1}, Float64}
    select :: Array{Bool, 1}
    # for gradient based inversion
    J      :: AbstractArray
    W      :: SparseMatrixCSC
    res    :: Vector
    nfixed :: Int
    ρ      
end

function Line(d::Array{Float64, 1} ;useML=false, σ=1.0, calcjacobian=false)
    if isa(σ, Array)
        @assert length(d) == length(σ)
    else
        σ = σ*ones(length(d))    
    end
    select = .!isnan.(d[:])
    res = copy(d[select])
    if calcjacobian
        J = getA(d)
        Wdiag = 1 ./σ[select]
    else    
        J, Wdiag = zeros(0), zeros(0)
    end 
    W = sparse(diagm(Wdiag))
    nfixed = 0 # only needed for gradient based
    Line(d, useML, σ, select, J, W, res, nfixed, copy(d))
end

function getA(v::AbstractVector)
    m = length(v)
    n = sum(.!isnan.(v))
    # @assert n<m
    sparse(1:n,findall(.!isnan.(v)),ones(n),n,m)
end  

function get_misfit(m::Model, opt::Options, line::Line)
    chi2by2 = 0.0
    if !opt.debug
        getr!(line, m.fstar[line.select])
        res, σ, select = line.res, line.σ, line.select     
        r = res./σ[select]
        chi2 = r'*r
        if line.useML
            N = sum(select)
            chi2by2 = 0.5*N*log(chi2)
        else
            chi2by2 = 0.5*chi2
        end
    end
    return chi2by2
end

function getr!(line, m)
    d, σ, select = line.d, line.σ, line.select
    line.res[:] = (m - d[select])
    nothing
end

function getresidual(line::Line, m::Vector{Float64}; computeJ=false)
    getr!(line, m[line.select])
end

function makehist(line::Line, opt::Options;
    nbins=100, burninfrac=0.5, temperaturenum=1)
    linidx = findall(.!isnan.(line.d))
    M = assembleTat1(opt, :fstar,
        temperaturenum=temperaturenum, burninfrac=burninfrac)
    resids = similar(M)
    for (im, m) in enumerate(M)
        resids[im] = m[sort(linidx)] - line.d[sort(linidx)]
    end
    resids = permutedims(hcat(resids...))
    mmin, mmax = extrema(resids)
    edges = LinRange(mmin, mmax, nbins+1)
    himage = zeros(Float64, nbins,length(linidx))
    for ilayer=1:length(linidx)
        himage[:,ilayer] = fit(Histogram, resids[:,ilayer], edges).weights
        himage[:,ilayer] = himage[:,ilayer]/sum(himage[:,ilayer])/(diff(edges)[1])
    end
    himage, edges
end

end
