module GP
using Statistics, LinearAlgebra, Distributions
function gaussiankernel(x::AbstractArray, y::AbstractArray, p)
    return exp(-0.5*norm(x-y,p)^p)
end

function makekernel(xtrain::AbstractArray, xtest::AbstractArray, λ::AbstractArray, p)
    # each input data arranged in a column
    nrows = size(xtrain, 2)
    ncols = size(xtest, 2)
    ϕ = zeros(nrows, ncols)
    tempx, tempy = zeros(Float64, size(xtrain,1)), zeros(Float64, size(xtrain,1))
    for i = 1:nrows, j = 1:ncols
        for k = 1:size(xtrain,1)
            tempx[k] = xtrain[k,i]/λ[k]
            tempy[k] =  xtest[k,j]/λ[k]
        end
            ϕ[i,j] = gaussiankernel(tempx, tempy, p)
    end
    return ϕ
end

function GPfit(ytrain, xtrain, xtest, λ::Array{Float64,1}, δ ;nogetvars=false, demean=true, p=2)
    @assert length(λ) == size(xtrain,1)
    my = 0
    if demean
        my = mean(ytrain)
    end
    ytrain = ytrain .- my
    Kstar = GP.makekernel(xtrain, xtest, λ, p)
    K_y = GP.makekernel(xtrain, xtrain, λ, p) + Matrix(δ^2*I, (size(xtrain,2)), (size(xtrain,2)))
    U = cholesky(K_y).U
    ytest = my .+ Kstar'*(U\(U'\ytrain))
    var_prior, var_y = [],[]
    if !nogetvars
        var_prior = GP.makekernel(xtest, xtest, λ, p)
        var_y =  var_prior - Kstar'*(U\(U'\Kstar))
    end
    ytest, var_y, var_prior
end

function makekernel(xtrain::AbstractArray, xtest::AbstractArray, λtest::Array{Float64, 2},
    λtrain::Array{Float64,2}, p)
    # each input data arranged in a column
    nrows = size(xtrain, 2)
    ncols = size(xtest, 2)
    ϕ = zeros(nrows, ncols)
    tempx, tempy = zeros(Float64, size(xtrain,1)), zeros(Float64, size(xtrain,1))
    for i = 1:nrows, j = 1:ncols
        # for k = 1:size(xtrain,1)
        #     tempx[k] = xtrain[k,i]/λ[k]
        #     tempy[k] =  xtest[k,j]/λ[k]
        # end
        #     ϕ[i,j] = gaussiankernel(tempx, tempy, p)
        ϕ[i,j] = 2^(size(λtrain,1)/2)*
                (det(diagm(0=>λtrain[:,i].^2))*det(diagm(0=>λtest[:,j].^2)))^0.25 *
                det(diagm(0=>(λtrain[:,i].^2 + λtest[:,j].^2)))^-0.5 *exp(-0.5*((sqrt(2)*norm((xtrain[:,i] - xtest[:,j])./(sqrt.(λtrain[:,i].^2 + λtest[:,j].^2)),p))^p))
    end
    return ϕ
end

function GPfit(ytrain, xtrain, xtest, λtest::Array{Float64,2}, λtrain::Array{Float64,2},
            δ ;nogetvars=false, demean=true, p=2)
    @assert size(λtest) == size(xtest)
    @assert size(λtrain) == size(xtrain)
    my = 0
    if demean
        my = mean(ytrain)
    end
    ytrain = ytrain .- my
    Kstar = GP.makekernel(xtrain, xtest, λtest, λtrain, p)
    K_y = GP.makekernel(xtrain, xtrain, λtrain, λtrain, p) + Matrix(δ^2*I, (size(xtrain,2)), (size(xtrain,2)))
    U = cholesky(K_y).U
    ytest = my .+ Kstar'*(U\(U'\ytrain))
    var_prior, var_y = [],[]
    if !nogetvars
        var_prior = GP.makekernel(xtest, xtest, λtest, p)
        var_y =  var_prior - Kstar'*(U\(U'\Kstar))
    end
    ytest, var_y, var_prior
end

function getPI(y::Array{Float64, 1}, fₓ::Array{Float64, 1}, σₓ²::Array{Float64, 1}; tol= 1e-12)
    d = Normal()
    X = (maximum(y) .- fₓ .- tol)./sqrt.(σₓ²)
    return ccdf.(d, X)
end

end
