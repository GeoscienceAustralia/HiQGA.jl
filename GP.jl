module GP
using Statistics, LinearAlgebra
function gaussiankernel(x::AbstractArray, y::AbstractArray, p)
    return exp(-0.5*norm(x-y,p)^p)
end

abstract type Kernel end

struct SqEuclidean <: Kernel end
struct Mat32 <:Kernel end
struct Mat52 <: Kernel end

κ(K::SqEuclidean, d::Real, p::Real) = exp(-0.5*(d^p))
κ(K::Mat32, d::Real, p::Real) = (1.0 + sqrt(3.0) * d) * exp(-sqrt(3.0) * d)
κ(K::Mat52, p::Real) = (1.0 + sqrt(5.0) * d + 5.0 * d^2 / 3.0) * exp(-sqrt(5.0) * d)

function makekernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                    λ::AbstractArray, p::Real)
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
            d = norm(tempx-tempy, p)
            ϕ[i,j] = κ(K, d, p)
    end
    return ϕ
end

function GPfit(K::Kernel, ytrain, xtrain, xtest, λ::Array{Float64,1}, δ::Real ;nogetvars=false,
            demean=true, p=2)
    @assert length(λ) == size(xtrain,1)
    my = 0
    if demean
        my = mean(ytrain)
    end
    ytrain = ytrain .- my
    Kstar = GP.makekernel(K, xtrain, xtest, λ, p)
    K_y = GP.makekernel(K, xtrain, xtrain, λ, p) + Matrix(δ^2*I, (size(xtrain,2)), (size(xtrain,2)))
    U = cholesky(K_y).U
    ytest = my .+ Kstar'*(U\(U'\ytrain))
    var_prior, var_y = [],[]
    if !nogetvars
        var_prior = GP.makekernel(K, xtest, xtest, λ, p)
        var_y =  var_prior - Kstar'*(U\(U'\Kstar))
    end
    ytest, var_y, var_prior
end

function makekernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
    λtest::Array{Float64, 2}, λtrain::Array{Float64,2}, p)
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
        d = norm((xtrain[:,i] - xtest[:,j])./(sqrt.(0.5*(λtrain[:,i].^2 + λtest[:,j].^2))), p)
        ϕ[i,j] = 2^(size(λtrain,1)/2)*
                (det(diagm(0=>λtrain[:,i].^2))*det(diagm(0=>λtest[:,j].^2)))^0.25 *
                det(diagm(0=>(λtrain[:,i].^2 + λtest[:,j].^2)))^-0.5 *
                κ(K, d, p)
    end
    return ϕ
end

function GPfit(K::Kernel, ytrain, xtrain, xtest, λtest::Array{Float64,2}, λtrain::Array{Float64,2},
            δ ;nogetvars=false, demean=true, p=2)
    @assert size(λtest) == size(xtest)
    @assert size(λtrain) == size(xtrain)
    my = 0
    if demean
        my = mean(ytrain)
    end
    ytrain = ytrain .- my
    Kstar = GP.meshkernel(K, xtrain, xtest, λtest, λtrain, p)
    K_y = GP.meshkernel(K, xtrain, xtrain, λtrain, λtrain, p) + Matrix(δ^2*I, (size(xtrain,2)), (size(xtrain,2)))
    U = cholesky(K_y).U
    ytest = my .+ Kstar*(U\(U'\ytrain))
    var_prior, var_y = [],[]
    if !nogetvars
        var_prior = GP.meshkernel(K, xtest, xtest, λtest, λtest, p)
        var_y =  var_prior - Kstar*(U\(U'\Kstar'))
    end
    ytest, var_y, var_prior
end

function getPI(y::Array{Float64, 1}, fₓ::Array{Float64, 1}, σₓ²::Array{Float64, 1}; tol= 1e-12)
    d = Normal()
    X = (maximum(y) .- fₓ .- tol)./sqrt.(σₓ²)
    return ccdf.(d, X)
end

function meshkernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                    λtest::AbstractArray, λtrain::AbstractArray, p)
    nrows = size(xtest, 2)
    ncols = size(xtrain, 2)
        ϕ = [kernel(K, xtrain[:,j], xtest[:,i], λtrain[:,j], λtest[:,i]) for i = 1:nrows, j = 1:ncols]
    return ϕ
end

function mapkernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                    λtest::AbstractArray, λtrain::AbstractArray; p=2)
    nrows = size(xtest, 2)
    ncols = size(xtrain, 2)
    ϕ = Array{Float64, 2}(undef, nrows, ncols)
    map!(x->x, ϕ, pairwise(K, xtrain, xtest, λtrain, λtest))
end

function colwise(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                 λtrain::AbstractArray, λtest::AbstractArray)
    na, nb = size(xtrain, 2), size(xtest, 2)
    @assert min(na, nb) == 1
    if na == 1
        a, b, la, lb = xtrain, xtest, λtrain, λtest
    else
        a, b, la, lb = xtest, xtrain, λtest, λtrain
    end
    [kernel(K::Kernel, a, b[:,j], la, lb[:,j]) for j = 1:max(na, nb)]
end

function pairwise(K::Kernel, xtrain, xtest, λtrain, λtest)
    nrows = size(xtest, 2)
    ncols = size(xtrain, 2)
    [kernel(K, xtrain[:,j], xtest[:,i], λtrain[:,j], λtest[:,i]) for i = 1:nrows, j = 1:ncols]
end

function kernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                λtest::AbstractArray, λtrain::AbstractArray; p=2)
    avλ² = 0.5*(λtrain.^2 + λtest.^2)
    dist = norm((xtrain - xtest)./sqrt.(avλ²),p)
    det(diagm(0=>(λtrain.*λtest).^2))^0.25 *
    det(diagm(0=>(avλ²)))^-0.5 *
    κ(K, dist, p)
end

end
