using LinearMaps, SparseArrays, PositiveFactorizations
using .AbstractOperator, .GP

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

function makereg(r::Symbol, F::Operator)
    r == :R0 && return sparse(makeregR0(F))
    r == :R1 && return sparse(makeregR1(F))
    r == :R2 && return sparse(makeregR1(F)*makeregR1(F))
    error("unknown regularization")
end   

function pushback(m, lo, hi)
    for i in eachindex(m)
        while (m[i]<lo) || (m[i]>hi)
                (m[i]<lo) && (m[i] = 2*lo - m[i])
                (m[i]>hi) && (m[i] = 2*hi - m[i])
        end
    end
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
                   F::Operator, λ²::Vector{Float64}, R::SparseMatrixCSC, target, lo, hi;
                   regularizeupdate = false)
    getresidual(F, m, computeJ=true)
    r, W = F.res, F.W
    @info "χ² is $(norm(W*r)^2)"
    for (i, l²) in enumerate(λ²)
        mnew[i] = m + newtonstep(m, m0, F, l², R, 
                        regularizeupdate=regularizeupdate)
        pushback(mnew[i], lo, hi)                
        getresidual(F, mnew[i], computeJ=false)
        push!(χ², norm(W*r)^2)
    end
    idx = -1
    if all(χ² .> target)
        idx = argmin(χ²)
    else    
        idx = findlast(χ² .<= target)
    end
    idx
end 

function bostep(m::AbstractVector, m0::AbstractVector, mnew::Vector{Vector{Float64}}, χ²::Vector{Float64},
                   F::Operator, λ²::Vector{Float64}, R::SparseMatrixCSC, target, lo, hi;
                   regularizeupdate = false, 
                   ## GP stuff
                   demean = true, 
                   κ = GP.SqEuclidean(), 
                   λ²GP = NaN, 
                   δtry = 1e-2,
                   frac = 5,
                   ntestdivs = 50,
                   acqfun = GP.EI(),
                   ntries = 6,
                   knownvalue = 0.7)
    
    getresidual(F, m, computeJ=true)
    r, W = F.res, F.W
    χ²₀ = norm(W*r)^2
    @info "χ² is $χ²₀ at start"
    
    knownvalue *= χ²₀

    lmin, lmax = extrema(log10.(λ²))
    λ²GP = [((lmax-lmin)/frac)^2] # length scale square for surrogate
    t = LinRange(lmin, lmax, ntestdivs) # test range for surrogate
    ttrain = [t[1]] # first training point location
    mnew[1] = m + newtonstep(m, m0, F, 10^ttrain[end], R, 
                        regularizeupdate=regularizeupdate)
    pushback(mnew[1], lo, hi)
    getresidual(F, mnew[1], computeJ=false)                
    @info "χ² is $(norm(W*r)^2) after first push"
    push!(χ², norm(W*r)^2) # first training value
    ytest, σ2, = GP.GPfit(κ, χ²', ttrain', t', λ²GP, δtry, demean=demean) # fit a GP
    # @info χ², ttrain
    # @info vec(ytest)
    # @info diag(σ2)
    AF = GP.getAF(acqfun, χ², vec(ytest), diag(σ2), findmin=true, knownvalue=knownvalue) # acquisition func calculation
    # @info AF
    # @info "---"
    for i = 2:ntries
        nextpos = argmax(AF) # get next test location
        # find the misfit
        mnew[i] = m + newtonstep(m, m0, F, 10^t[nextpos], R, 
        regularizeupdate=regularizeupdate)
        pushback(mnew[i], lo, hi)                        
        getresidual(F, mnew[i], computeJ=false)
        @show push!(χ², norm(W*r)^2) # next training value
        push!(ttrain, t[nextpos]) # next training location
        ytest, σ2, = GP.GPfit(κ, χ²', ttrain', t', λ²GP, δtry, demean=demean)
        # @info χ², ttrain
        # @info vec(ytest)
        # @info diag(σ2)
        AF = GP.getAF(acqfun, χ², vec(ytest), diag(σ2), findmin=true, knownvalue=knownvalue)
        # @info AF
    end    
    idx = -1 
    if all(χ² .> target)
        idx = argmin(χ²)
    else    
        idx = findlast(χ² .<= target)
    end
    idx
end 

function gradientinv(   m::AbstractVector,
                        m0::AbstractVector, 
                        F::Operator, λ²::Vector{Float64}; 
                        regtype=:R0,
                        saveall = true, 
                        nstepsmax = 10,
                        target = nothing,
                        lo = -3.,
                        hi = 1.,
                        regularizeupdate = false,
                        dobo=true)
    R = makereg(regtype, F)                
    ndata = length(F.res)
    isnothing(target) && (target = ndata)
    if saveall
        mnew = [[similar(m) for i in 1:length(λ²)] for j in 1:nstepsmax]
        χ²   = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
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
        if dobo
            idx = bostep(m, m0, mn, χsq, F, λ², R, ndata, lo, hi,
            regularizeupdate=regularizeupdate)
        else                
            idx = occamstep(m, m0, mn, χsq, F, λ², R, ndata, lo, hi,
                        regularizeupdate=regularizeupdate)
        end                
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