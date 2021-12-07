module GP

using Statistics, LinearAlgebra, Distributions
function gaussiankernel(x::AbstractArray, y::AbstractArray, p)
    return exp(-0.5*norm(x-y,p)^p)
end

abstract type Kernel end

struct SqEuclidean <: Kernel end
struct Mat32 <:Kernel end
struct Mat52 <: Kernel end
struct OrstUhn <: Kernel end

κ(K::SqEuclidean, d::Real; p::Real=2.0) = exp(-0.5*(d^p))
κ(K::Mat32, d::Real; p::Real=2.0) = (1.0 + sqrt(3.0) * d) * exp(-sqrt(3.0) * d)
κ(K::Mat52, d::Real; p::Real=2.0 ) = (1.0 + sqrt(5.0) * d + 5.0 * abs2(d) / 3.0) * exp(-sqrt(5.0) * d)
κ(K::OrstUhn, d::Real; p::Real=2.0) = exp(-d)

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
            demean=true, p=2, my = zeros(size(ytrain, 1)))
    @assert length(λ²) == size(xtrain,1)
    if demean
        my = mean(ytrain, dims=2)
    end
    ytrain = ytrain .- my
    Kstar = makekernel(K, xtrain, xtest, λ², p)
    K_y = makekernel(K, xtrain, xtrain, λ², p) + Matrix(δ^2*I, (size(xtrain,2)), (size(xtrain,2)))
    U = cholesky(K_y).U
    ytest = my' .+ Kstar'*(U\(U'\ytrain'))
    var_prior, var_y = [],[]
    if !nogetvars
        var_prior = makekernel(K, xtest, xtest, λ², p)
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
            δ ;nogetvars=false, demean=true, p=2, my = zeros(size(ytrain, 1)))
    @assert size(λ²test) == size(xtest)
    @assert size(λ²train) == size(xtrain)
    if demean
        my = mean(ytrain, dims=2)
    end
    ytrain = ytrain .- my
    Kstar = meshkernel(K, xtrain, xtest, λ²test, λ²train, p)
    K_y = meshkernel(K, xtrain, xtrain, λ²train, λ²train, p) + Matrix(δ^2*I, (size(xtrain,2)), (size(xtrain,2)))
    U = cholesky(K_y).U
    ytest = my' .+ Kstar*(U\(U'\ytrain'))
    var_prior, var_y = [],[]
    if !nogetvars
        var_prior = meshkernel(K, xtest, xtest, λ²test, λ²test, p)
        var_y =  var_prior - Kstar*(U\(U'\Kstar'))
    end
    ytest, var_y, var_prior
end

abstract type AcqFunc end
struct EI <: AcqFunc end
struct PI <: AcqFunc end

function getAF(AF::PI, y::Array{Float64, 1}, fₓ::Array{Float64, 1}, σₓ²::Array{Float64, 1}, 
        tol = 1e-12, kv = NaN)
    isnan(kv) && (kv = maximum(y))
    d = Normal()
    X = (fₓ .- kv .- tol)./sqrt.(σₓ²)
    cdf.(d, X)
end
    
function getAF(AF::EI, y::Array{Float64, 1}, fₓ::Array{Float64, 1}, σₓ²::Array{Float64, 1}, 
        tol = 1e-12, kv = NaN)
    isnan(kv) && (kv = maximum(y))    
    d = Normal()
    X = (fₓ .- kv .- tol)
    σₓ = sqrt.(σₓ²)
    Z = X./σₓ
    X.*cdf.(d, Z) + σₓ.*pdf.(d, Z) 
end

function getAF(AF::AcqFunc, y::Vector{Float64}, fₓ::Vector{Float64}, σₓ²::Vector{Float64};
    findmin = false, tol = 1e-12, knownvalue=NaN)
    sf = length(y) == 1 ? 1 : var(y) 
    if findmin
        AF = getAF(AF, -y, -fₓ, sf*σₓ², tol, -knownvalue)
    else
        AF = getAF(AF,  y,  fₓ, sf*σₓ², tol, knownvalue)    
    end
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
    # @views begin
        [@inbounds @fastmath kernel(K::Kernel, a, b[:,j], la, lb[:,j]) for j = 1:max(na, nb)]
    # end
end

function colwise!(A::AbstractArray, K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                 λ²train::AbstractArray, λ²test::AbstractArray)
    na, nb = size(xtrain, 2), size(xtest, 2)
    @assert min(na, nb) == 1
    if na == 1
        a, b, la, lb = xtrain, xtest, λ²train, λ²test
    else
        a, b, la, lb = xtest, xtrain, λ²test, λ²train
    end
    @views begin
        @inbounds @fastmath for j = 1:max(na, nb)
            A[j] = kernel(K::Kernel, a, b[:,j], la, lb[:,j])
        end
    end
end

function pairwise(A::AbstractArray, K::Kernel, xtrain::AbstractArray, xtest::AbstractArray, λ²train::AbstractArray, λ²test::AbstractArray)
    nrows = size(xtest, 2)
    ncols = size(xtrain, 2)
    @views begin
        @inbounds @fastmath for j = 1:ncols, i = 1:nrows
            A[i,j] = kernel(K, xtrain[:,j], xtest[:,i], λ²train[:,j], λ²test[:,i])
        end
    end
end

@inline function kernel(K::Kernel, xtrain::AbstractArray, xtest::AbstractArray,
                λ²test::AbstractArray, λ²train::AbstractArray; p=2)
    t = eltype(xtrain)
    d, c = zero(t), one(t)
    @simd for i = 1:length(xtrain)
        avλ² = 0.5*(λ²train[i] + λ²test[i])
        d += abs2(xtrain[i] - xtest[i])/avλ²
        c *= sqrt(sqrt((λ²train[i]*λ²test[i])/(avλ²*avλ²)))
    end
    c*κ(K, sqrt(d), p=p)
end

end
