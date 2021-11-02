using LinearMaps, SparseArrays, PositiveFactorizations, Roots
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

function occamstep(m::AbstractVector, m0::AbstractVector, mnew::Vector{Vector{Float64}},
                χ²::Vector{Float64}, λ²::Vector{Vector{Float64}}, F::Operator, R::SparseMatrixCSC, target, 
                lo, hi, λ²min, λ²max, ntries; knownvalue=NaN, regularizeupdate = false, ϵ=.001)
                   
    getresidual(F, m, computeJ=true)
    r, W = F.res, F.W
    r₀ = copy(r)    
    χ²₀ = norm(W*r)^2
    knownvalue *= χ²₀
    α = 1
    count = 0
    countmax = 6
    while true
        for (i, l²) in enumerate(10 .^LinRange(λ²min, λ²max, ntries))
            if count == 0    
                mnew[i] = m + newtonstep(m, m0, F, l², R, 
                                regularizeupdate=regularizeupdate)
                pushback(mnew[i], lo, hi)                
                getresidual(F, mnew[i], computeJ=false)
                push!(χ², norm(W*r)^2)
                push!(λ², [l²; α])
                r .= r₀
            else    
                mnew[i] = α*mnew[i] + (1 - α)*m
                pushback(mnew[i], lo, hi)                
                getresidual(F, mnew[i], computeJ=false)
                χ²[i] = norm(W*r)^2
                λ²[i] = [l²; α]
                r .= r₀
            end
            χ²[i] <= knownvalue && (count = countmax; break)
        end
        if all(χ² .>= χ²₀)
            α = α/2
        else 
            count = countmax    
        end
        count += 1        
        count > countmax && break                    
    end
    idx = -1
    if all(χ² .> target)
        idx = argmin(χ²)
    else    
        idx = findlast(χ² .<= target)
        if 1 - χ²[idx]/target >  - ϵ
            # # root finding
            # @info "before $(λ²[idx][1])"
            # λ²[idx][1] = find_zero(l² -> fλ²(α, m, m0, F, l², R, regularizeupdate, lo, hi) - target, (λ²[idx][1], λ²[idx+1][1]))
            # @info "after $(λ²[idx][1])"
            # mnew[idx] = m + α*newtonstep(m, m0, F, λ²[idx][1], R, regularizeupdate=regularizeupdate)
            # χ²[idx] = getχ²(F, mnew[idx])
            dophase2(m, m0, F, R, regularizeupdate, lo, hi, target, idx, idx+1, mnew,  χ², λ²)
        end    
    end
    idx
end 

function fλ²(α, m, m0, F, l², R, regularizeupdate, lo, hi)
    mnew = m + α*newtonstep(m, m0, F, l², R, regularizeupdate=regularizeupdate)
    pushback(mnew, lo, hi)                
    getχ²(F, mnew)
end

function getχ²(F, m)
    r, W = F.res, F.W
    r₀ = copy(r)
    getresidual(F, m, computeJ=false)
    χ² = norm(W*r)^2
    r .= r₀
    return χ²
end

function initbo(λ²min, λ²max, λ²frac, ntestdivsλ², αmin, αmax, αfrac, ntestdivsα)
    l2 = LinRange(λ²min, λ²max, ntestdivsλ²) # test range for surrogate
    α = LinRange(αmin, αmax, ntestdivsα) # test range for surrogate
    λ²GP = [((λ²max-λ²min)/λ²frac)^2; (abs(α[end]-α[1])/αfrac)^2] # length scale square for surrogate
    X1, X2 = [x1 for x1 in l2, x2 in α], [x2 for x1 in l2, x2 in α]
    t = [X1[:]';X2[:]']
    t, λ²GP
end    

function bostep(m::AbstractVector, m0::AbstractVector, mnew::Vector{Vector{Float64}}, χ²::Vector{Float64}, λ²sampled::Vector{Vector{Float64}},
                    t::Array{Float64, 2}, λ²GP::Array{Float64, 1}, F::Operator, R::SparseMatrixCSC, target, lo, hi;
                   regularizeupdate = false, ϵ=.001,
                  
                   ## GP stuff
                   demean = true, 
                   κ = GP.Mat52(), 
                   δtry = 1e-2,
                   acqfun = GP.EI(),
                   ntries = 6,
                   firstvalue = :last,
                   knownvalue = NaN,
                   breakonknown = false)
    
    getresidual(F, m, computeJ=true)
    r, W = F.res, F.W
    r₀ = copy(r)    
    χ²₀ = norm(W*r)^2
    knownvalue *= χ²₀

    ttrain = zeros(size(t,1), 0)
    for i = 1:ntries
        nextpos, = getBOsample(κ, χ²', ttrain, t, λ²GP, δtry, demean, i, knownvalue, firstvalue, acqfun)
        push!(λ²sampled, [10^t[1, nextpos]; 2^t[2,nextpos]])
        mnew[i] = m + 2^t[2,nextpos]*newtonstep(m, m0, F, 10^t[1,nextpos], R, regularizeupdate=regularizeupdate)
        pushback(mnew[i], lo, hi)                        
        getresidual(F, mnew[i], computeJ=false)
        push!(χ², norm(W*r)^2) # next training value
        r .= r₀
        ttrain = hcat(ttrain, t[:,nextpos]) # next training location
        (χ²[i] <= knownvalue && breakonknown) && break
    end    
    idx = -1 
    if all(χ² .> target)
        idx = argmin(χ²)
    else 
        sortedλ²idx = sortperm(vec(ttrain[1,:]))   
        sortedχ²idx = findlast(χ²[sortedλ²idx] .<= target)
        idx = sortedλ²idx[sortedχ²idx]
        if 1 - χ²[idx]/target >  - ϵ
            # root finding
            # closing bracket guess index
            b = sortedλ²idx[sortedχ²idx+1]
            dophase2(m, m0, F, R, regularizeupdate, lo, hi, target, idx, b, mnew,  χ², λ²sampled)
        end
    end
    idx
end 

function dophase2(m, m0, F, R, regularizeupdate, lo, hi, target, aidx, bidx, mnew,  χ², λ²)
    α = λ²[aidx][2]
    λ²[aidx][1] = find_zero(l² -> fλ²(α, m, m0, F, l², R, regularizeupdate, lo, hi) - target, (λ²[aidx][1], λ²[bidx][1]))
    mnew[aidx] = m + α*newtonstep(m, m0, F, λ²[aidx][1], R, regularizeupdate=regularizeupdate)
    χ²[aidx] = getχ²(F, mnew[aidx])
    nothing
end    

function getBOsample(κ, χ², ttrain, t, λ²GP, δtry, demean, iteration, knownvalue, firstvalue, acqfun)
    # χ², ttrain, t are row major
    ntrain = length(ttrain)
    if ntrain > 0
        ytest, σ2, = GP.GPfit(κ, χ², ttrain, t, λ²GP, δtry, demean=demean)
        diagσ2 = diag(σ2)
        nextpos = argmax(GP.getAF(acqfun, vec(χ²), vec(ytest), diagσ2, findmin=true, knownvalue=knownvalue))          
    else
        nextpos = AFoneiter(firstvalue, size(t, 2))
    end
    nextpos
end    

function AFoneiter(firstvalue::Symbol, n::Int)
    if firstvalue == :first
        return 1
    elseif firstvalue == :last
        return n
    elseif firstvalue == :middle
        return round(Int, middle(1:n))
    elseif firstvalue == :random
        return rand(1:n)
    end
end

function gradientinv(   m::AbstractVector,
                        m0::AbstractVector, 
                        F::Operator; 
                        regtype=:R0,
                        nstepsmax = 10,
                        ntries = 6,
                        target = nothing,
                        lo = -3.,
                        hi = 1.,
                        λ²min = 0,
                        λ²max = 8,
                        λ²frac=5,
                        ntestdivsλ²=50,
                        αmin=-4, 
                        αmax=0, 
                        αfrac=4, 
                        ntestdivsα=32,
                        regularizeupdate = false,
                        knownvalue=0.7,
                        firstvalue=:last,
                        κ = GP.Mat52(),
                        breakonknown=false,
                        dobo = false)
    R = makereg(regtype, F)                
    ndata = length(F.res)
    isnothing(target) && (target = ndata)
    mnew = [[similar(m) for i in 1:ntries] for j in 1:nstepsmax]
    χ²   = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
    λ² = [Vector{Vector{Float64}}(undef, 0) for j in 1:nstepsmax]
    oidx = zeros(Int, nstepsmax)  
    ndata = length(F.res)
    t, λ²GP = initbo(λ²min, λ²max, λ²frac, ntestdivsλ², αmin, αmax, αfrac, ntestdivsα)
    istep = 1                  
    while true
        if dobo
            idx = bostep(m, m0, mnew[istep], χ²[istep], λ²[istep], t, λ²GP, F, R, ndata, lo, hi,
            regularizeupdate=regularizeupdate, ntries=ntries, κ = κ,
            knownvalue=knownvalue, firstvalue=firstvalue, breakonknown=breakonknown)         
        else
            idx = occamstep(m, m0, mnew[istep], χ²[istep], λ²[istep], F, R, target, 
                lo, hi, λ²min, λ²max, ntries, knownvalue=knownvalue, regularizeupdate = regularizeupdate)
        end
        @info "iteration: $istep χ²: $(χ²[istep][idx]) target: $target"
        m = mnew[istep][idx]
        oidx[istep] = idx
        χ²[istep][idx] < target && break
        istep += 1
        istep > nstepsmax && break
    end
    return mnew, χ², λ², oidx
end    