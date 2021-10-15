using LinearMaps, SparseArrays
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

function newtonstep(m::AbstractVector, m0::AbstractVector, F::Operator, λ²::Float64, R::LinearMap; 
                    regularizeupdate=true)
    J, r, Cdinv = F.J, F.res, F.Cdinv
    if !regularizeupdate # regularize model
        -sparse(J'Cdinv*J + λ²*R'R)\(J'Cdinv*r + λ²*R'R*(m - m0))
    else
        -sparse(J'Cdinv*J + λ²*R'R)\(J'Cdinv*r)
    end       
end

function occamstep(m::AbstractVector, m0::AbstractVector, mnew::Vector{Vector{Float64}}, χ²::Vector{Float64},
                   F::Operator, λ²::Vector{Float64}, R::LinearMap, ndata;
                   regularizeupdate = false)
    getresidual(F, m, computeJ=true)
    r = F.res
    @info "χ² is $(r'F.Cdinv*r)"
    for (i, l²) in enumerate(λ²)
        mnew[i] = m + newtonstep(m, m0, F, l², R, 
                        regularizeupdate=regularizeupdate)
        getresidual(F, mnew[i], computeJ=false)
        χ²[i] = r'F.Cdinv*r
    end
    idx = -1
    if all(χ² .> ndata)
        idx = argmin(χ²)
    else    
        idx = findlast(χ² .<= ndata)
    end
    idx
end

function gradientinv(   m::AbstractVector,
                        m0::AbstractVector, 
                        F::Operator, λ²::Vector{Float64}; 
                        R=makeregR0(F),
                        saveall = true, 
                        nstepsmax = 10,
                        regularizeupdate = false)
                        
    ndata = length(F.res)
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
        @info "iteration: $istep χ²: $(χsq[idx]) target: $ndata"
        m = mn[idx]
        saveall && (oidx[istep] = idx)
        istep += 1
        χsq[idx] < ndata && break
        istep > nstepsmax && break
    end
    if saveall
        return mnew, χ², oidx
    end      
end    