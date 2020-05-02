module GP
using Statistics, LinearAlgebra
function gaussiankernel(x::AbstractArray, y::AbstractArray, p)
    return exp(-0.5*norm(x-y,p)^p)
end

abstract type Kernel end

struct SqEuclidean <: Kernel end
struct Mat32 <:Kernel end
struct Mat52 <: Kernel end

κ(K::SqEuclidean, d::Real; p::Real=2.0) = exp(-0.5*(d^p))
κ(K::Mat32, d::Real; p::Real=2.0) = (1.0 + sqrt(3.0) * d) * exp(-sqrt(3.0) * d)
κ(K::Mat52, d::Real; p::Real=2.0 ) = (1.0 + sqrt(5.0) * d + 5.0 * d^2 / 3.0) * exp(-sqrt(5.0) * d)

function makekernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                    λ²::AbstractArray, p::Real)
    # each input data arranged in a column
    nrows = size(xtrain, 2)
    ncols = size(xtest, 2)
    ϕ = zeros(nrows, ncols)
    tempx, tempy = zeros(Float64, size(xtrain,1)), zeros(Float64, size(xtrain,1))
    for i = 1:nrows, j = 1:ncols
        for k = 1:size(xtrain,1)
            tempx[k] = xtrain[k,i]/λ²[k]^0.5
            tempy[k] =  xtest[k,j]/λ²[k]^0.5
        end
            d = norm(tempx-tempy, p)
            ϕ[i,j] = κ(K, d, p=p)
    end
    return ϕ
end

function GPfit(K::Kernel, ytrain, xtrain, xtest, λ²::Array{Float64,1}, δ::Real ;nogetvars=false,
            demean=true, p=2)
    @assert length(λ²) == size(xtrain,1)
    my = zeros(size(ytrain, 1))
    if demean
        my = mean(ytrain, dims=2)
    end
    ytrain = ytrain .- my
    Kstar = GP.makekernel(K, xtrain, xtest, λ², p)
    K_y = GP.makekernel(K, xtrain, xtrain, λ², p) + Matrix(δ^2*I, (size(xtrain,2)), (size(xtrain,2)))
    U = cholesky(K_y).U
    ytest = my' .+ Kstar'*(U\(U'\ytrain'))
    var_prior, var_y = [],[]
    if !nogetvars
        var_prior = GP.makekernel(K, xtest, xtest, λ², p)
        var_y =  var_prior - Kstar'*(U\(U'\Kstar))
    end
    ytest, var_y, var_prior
end

function makekernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
    λ²test::Array{Float64, 2}, λ²train::Array{Float64,2}, p)
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
        d = norm((xtrain[:,i] - xtest[:,j])./(sqrt.(0.5*(λ²train[:,i] + λ²test[:,j]))), p)
        ϕ[i,j] = 2^(size(λ²train,1)/2)*
                (det(diagm(0=>λ²train[:,i]))*det(diagm(0=>λ²test[:,j])))^0.25 *
                det(diagm(0=>(λ²train[:,i] + λ²test[:,j])))^-0.5 *
                κ(K, d, p=p)
    end
    return ϕ
end

function GPfit(K::Kernel, ytrain, xtrain, xtest, λ²test::Array{Float64,2}, λ²train::Array{Float64,2},
            δ ;nogetvars=false, demean=true, p=2)
    @assert size(λ²test) == size(xtest)
    @assert size(λ²train) == size(xtrain)
    my = zeros(size(ytrain, 1))
    if demean
        my = mean(ytrain, dims=2)
    end
    ytrain = ytrain .- my
    Kstar = GP.meshkernel(K, xtrain, xtest, λ²test, λ²train, p)
    K_y = GP.meshkernel(K, xtrain, xtrain, λ²train, λ²train, p) + Matrix(δ^2*I, (size(xtrain,2)), (size(xtrain,2)))
    U = cholesky(K_y).U
    ytest = my' .+ Kstar*(U\(U'\ytrain'))
    var_prior, var_y = [],[]
    if !nogetvars
        var_prior = GP.meshkernel(K, xtest, xtest, λ²test, λ²test, p)
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
                    λ²test::AbstractArray, λ²train::AbstractArray, p)
    nrows = size(xtest, 2)
    ncols = size(xtrain, 2)
        ϕ = [kernel(K, xtrain[:,j], xtest[:,i], λ²train[:,j], λ²test[:,i], p=p) for i = 1:nrows, j = 1:ncols]
    return ϕ
end

function mapkernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                    λ²test::AbstractArray, λ²train::AbstractArray; p=2)
    nrows = size(xtest, 2)
    ncols = size(xtrain, 2)
    ϕ = Array{Float64, 2}(undef, nrows, ncols)
    map!(x->x, ϕ, pairwise(K, xtrain, xtest, λ²train, λ²test))
end

function colwise(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                 λ²train::AbstractArray, λ²test::AbstractArray)
    na, nb = size(xtrain, 2), size(xtest, 2)
    @assert min(na, nb) == 1
    if na == 1
        a, b, la, lb = xtrain, xtest, λ²train, λ²test
    else
        a, b, la, lb = xtest, xtrain, λ²test, λ²train
    end
    [kernel(K::Kernel, a, b[:,j], la, lb[:,j]) for j = 1:max(na, nb)]
end

function pairwise(K::Kernel, xtrain, xtest, λ²train, λ²test)
    nrows = size(xtest, 2)
    ncols = size(xtrain, 2)
    [kernel(K, xtrain[:,j], xtest[:,i], λ²train[:,j], λ²test[:,i]) for i = 1:nrows, j = 1:ncols]
end

function kernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                λ²test::AbstractArray, λ²train::AbstractArray; p=2)
    avλ² = 0.5*(λ²train + λ²test)
    dist = norm((xtrain - xtest)./sqrt.(avλ²),p)
    det(diagm(0=>(λ²train.*λ²test)))^0.25 *
    det(diagm(0=>(avλ²)))^-0.5 *
    κ(K, dist, p=p)
end

end
