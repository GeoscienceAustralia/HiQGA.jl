using LinearMaps, SparseArrays, PositiveFactorizations
using .AbstractOperator

function makeregR0(F::Operator1D)
    n = length(F.ρ) - F.nfixed
    LinearMap(I, n)
end

function makeregR1(F::Operator1D)
    n = length(F.ρ) - F.nfixed
    LinearMap(R1Dop, Rt1Dop, n)
end

function R1Dop(x::Vector)
    vcat(0, diff(x))
end

function Rt1Dop(y::Vector)
    x = vcat(-diff(y),y[end])
    x[1] = -y[2]
    x
end

function newtonstep(m::AbstractVector, m0::AbstractVector, F::Operator, λ²::Float64, R::SparseMatrixCSC; 
                    regularizeupdate=true)
    JtW, Wr = F.J'*F.W, F.W*F.res
    H = (JtW*(JtW)' + λ²*R'R)
    U = cholesky(Positive, H, Val{false}).U 
    if !regularizeupdate # regularize model
        -U\(U'\(JtW*Wr + λ²*R'R*(m - m0)))
    else
        -U\(U'\(JtW*Wr))
    end       
end

function occamstep(m::AbstractVector, m0::AbstractVector, mnew::Vector{Vector{Float64}}, χ²::Vector{Float64},
                   F::Operator, λ²::Vector{Float64}, R::SparseMatrixCSC, target;
                   regularizeupdate = false)
    getresidual(F, m, computeJ=true)
    r, W = F.res, F.W
    @info "χ² is $(norm(W*r)^2)"
    for (i, l²) in enumerate(λ²)
        mnew[i] = m + newtonstep(m, m0, F, l², R, 
                        regularizeupdate=regularizeupdate)
        getresidual(F, mnew[i], computeJ=false)
        χ²[i] = norm(W*r)^2
    end
    idx = -1
    if all(χ² .> target)
        idx = argmin(χ²)
    else    
        idx = findlast(χ² .<= target)
    end
    idx
end

function makereg(r::Symbol, F::Operator)
    r == :R0 && return sparse(makeregR0(F))
    r == :R1 && return sparse(makeregR1(F))
    r == :R2 && return sparse(makeregR1(F)*makeregR1(F))
    error("unknown regularization")
end    

function gradientinv(   m::AbstractVector,
                        m0::AbstractVector, 
                        F::Operator, λ²::Vector{Float64}; 
                        regtype=:R0,
                        saveall = true, 
                        nstepsmax = 10,
                        target = nothing,
                        regularizeupdate = false)
    R = makereg(regtype, F)                
    ndata = length(F.res)
    isnothing(target) && (target = ndata)
    if saveall
        mnew = [[similar(m) for i in 1:length(λ²)] for j in 1:nstepsmax]
        χ²   = [similar(λ²) for j in 1:nstepsmax]
        oidx = zeros(Int, nstepsmax)
    else
        mnew = [similar(m) for i in 1:length(λ²)]
        χ²   = [similar(λ²) for i in 1:length(λ²)]
    end        
    ndata = length(F.res)
    istep = 1                  
    while true
        if saveall
            mn = mnew[istep]
            χsq = χ²[istep]
        else 
            mn = mnew
            χsq = χ²
        end        
        idx = occamstep(m, m0, mn, χsq, F, λ², R, ndata, 
                        regularizeupdate=regularizeupdate)
        @info "iteration: $istep χ²: $(χsq[idx]) target: $target"
        m = mn[idx]
        saveall && (oidx[istep] = idx)
        istep += 1
        χsq[idx] < target && break
        istep > nstepsmax && break
    end
    if saveall
        return mnew, χ², oidx
    end      
end    