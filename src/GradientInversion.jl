using LinearMaps, SparseArrays, PositiveFactorizations, LazyGrids
using Roots:find_zero
using .AbstractOperator, .GP, Optim

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

function pushback(m, lo, hi, reflect=true)
    for i in eachindex(m)
        while (m[i]<lo) || (m[i]>hi)
                (m[i]<lo) && (m[i] = reflect ? 2*lo - m[i] : lo)
                (m[i]>hi) && (m[i] = reflect ? 2*hi - m[i] : hi)
        end
    end
end  

function newtonstep(m::AbstractVector, m0::AbstractVector, F::Operator, λ²::Float64, β²::Float64, R::SparseMatrixCSC; 
                    regularizeupdate=true)
    JtW, Wr = F.J'*F.W, F.W*F.res
    H = (JtW*(JtW)' + λ²*R'R + λ²*β²*I)
    U = cholesky(Positive, H, Val{false}).U 
    if !regularizeupdate # regularize model
        -U\(U'\(JtW*Wr + λ²*R'R*m + λ²*β²*(m - m0)))
    else
        -U\(U'\(JtW*Wr))
    end       
end

function occamstep(m::AbstractVector, m0::AbstractVector, Δm::AbstractVector, mnew::Vector{Vector{Float64}},
                χ²::Vector{Float64}, λ²::Vector{Vector{Float64}}, F::Operator, R::SparseMatrixCSC, target, 
                lo, hi, λ²min, λ²max, β², ntries; knownvalue=NaN, regularizeupdate = false)
                   
    χ²₀ = getχ²(F, m, computeJ=true)
    knownvalue *= χ²₀
    α = 1
    count = 0
    countmax = 6
    while true
        for (i, l²) in enumerate(10 .^reverse(LinRange(λ²min, λ²max, ntries)))
            count == 0 && (Δm[i] = newtonstep(m, m0, F, l², β², R, regularizeupdate=regularizeupdate))
            mnew[i] = m + α*Δm[i]
            pushback(mnew[i], lo, hi)                
            chi2 = getχ²(F, mnew[i])
            nu = [l²; α]
            count == 0 && (push!(χ², chi2); push!(λ², nu))
            χ²[i] = chi2
            λ²[i] = nu
            if χ²[i] <= knownvalue
                count > 0 && map(x->deleteat!(x, i+1:ntries), (mnew, χ², λ²))
                count = countmax
                break
            end    
        end
        if all(χ² .>= χ²₀)
           α = count < countmax - 1 ? α/2 : 0. # make sure we don't start going uphill again
        end    
        count += 1        
        count > countmax && break                    
    end
    idx, foundroot = -1, false
    if all(χ² .> target)
        idx = argmin(χ²)
    else    
        idx = findlast(χ² .<= target)
        a = log10(λ²[idx][1])
        b = λ²max # isa(bidx, Nothing) ? λ²max : log10(λ²[bidx][1]) # maybe just set default λ²max
        dophase2(m, m0, F, R, regularizeupdate, lo, hi, target, a, b, mnew,  χ², λ², β², idx)
        foundroot = true    
    end
    idx, foundroot
end 


function dophase2(m, m0, F, R, regularizeupdate, lo, hi, target, a, b, mnew,  χ², λ², β², idx, reflect=true)
    α = λ²[idx][2]
    if (fλ²(α, m, m0, F, 10. ^b, β², R, regularizeupdate, lo, hi) <= target)
        λ²[idx][1] = 10. ^b
    else    
        λ²[idx][1] = 10^find_zero(l² -> fλ²(α, m, m0, F, 10^l², β², R, regularizeupdate, lo, hi) - target, (a, b))
    end    
    mnew[idx] = m + α*newtonstep(m, m0, F, λ²[idx][1], β², R, regularizeupdate=regularizeupdate)
    pushback(mnew[idx], lo, hi, reflect)
    χ²[idx] = getχ²(F, mnew[idx])
    nothing
end    

function fλ²(α, m, m0, F, l², β², R, regularizeupdate, lo, hi, reflect=true)
    mnew = m + α*newtonstep(m, m0, F, l², β², R, regularizeupdate=regularizeupdate)
    pushback(mnew, lo, hi, reflect)                
    getχ²(F, mnew)
end

function getχ²(F, m; computeJ=false)
    r, W = F.res, F.W
    r₀ = copy(r)
    getresidual(F, m, computeJ=computeJ)
    χ² = norm(W*r)^2
    computeJ || (r .= r₀)
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
                    t::Array{Float64, 2}, β²::Float64, λ²GP::Array{Float64, 1}, F::Operator, R::SparseMatrixCSC, target, lo, hi;
                   regularizeupdate = false, 
                  
                   ## GP stuff
                   demean = true, 
                   κ = GP.Mat52(), 
                   δtry = 1e-2,
                   acqfun = GP.EI(),
                   ntries = 6,
                   firstvalue = :last,
                   knownvalue = NaN,
                   breakonknown = false)
    
    χ²₀ = getχ²(F, m, computeJ=true)
    knownvalue *= χ²₀               

    ttrain = zeros(size(t,1), 0)
    for i = 1:ntries
        nextpos, = getBOsample(κ, χ²', ttrain, t, λ²GP, δtry, demean, knownvalue, firstvalue, acqfun)
        push!(λ²sampled, [10^t[1, nextpos]; 2^t[2,nextpos]])
        mnew[i] = m + 2^t[2,nextpos]*newtonstep(m, m0, F, 10^t[1,nextpos], β², R, regularizeupdate=regularizeupdate)
        pushback(mnew[i], lo, hi)                        
        push!(χ², getχ²(F, mnew[i])) # next training value
        ttrain = hcat(ttrain, t[:,nextpos]) # next training location
        (χ²[i] <= knownvalue && breakonknown) && break
    end
    idx, foundroot = -1, false 
    if all(χ² .> target)
        idx = argmin(χ²)
    else 
        sortedλ²idx = sortperm(vec(ttrain[1,:]))   
        sortedχ²idx = findlast(χ²[sortedλ²idx] .<= target)
        idx = sortedλ²idx[sortedχ²idx]
        a = log10(λ²sampled[idx][1])
        b = maximum(t[1,:]) # λ²max in log10 units
        dophase2(m, m0, F, R, regularizeupdate, lo, hi, target, a, b, mnew,  χ², λ²sampled, β², idx)
        foundroot = true
    end
    idx, foundroot
end 

function getBOsample(G:: GP.KernelStruct, χ², ttrain, t, knownvalue, firstvalue, acqfun)
    # χ², ttrain, t are row major
    ntrain = length(ttrain)
    if ntrain > 0
        # ytest, σ2, = GP.GPfit(κ, χ², ttrain, t, λ²GP, δtry, demean=demean)
        GP.GPfitaddpoint(G, ttrain, χ², t)
        ytest, σ2 = G.ytest, G.var_post
        nextpos = argmax(GP.getAF(acqfun, vec(χ²), vec(ytest), diag(σ2), findmin=true, knownvalue=knownvalue))          
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
                        regtype=:R1,
                        nstepsmax = 10,
                        ntries = 6,
                        target = nothing,
                        lo = -3.,
                        hi = 1.,
                        λ²min = 0,
                        λ²max = 8,
                        λ²frac=5,
                        β² = 0.,
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
                        dobo = false,
                        fname="")
    R = makereg(regtype, F)                
    ndata = length(F.res)
    isnothing(target) && (target = ndata)
    target₀ = target
    mnew = [[similar(m) for i in 1:ntries] for j in 1:nstepsmax]
    χ²   = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
    λ² = [Vector{Vector{Float64}}(undef, 0) for j in 1:nstepsmax]
    oidx = zeros(Int, nstepsmax)  
    ndata = length(F.res)
    t, λ²GP = initbo(λ²min, λ²max, λ²frac, ntestdivsλ², αmin, αmax, αfrac, ntestdivsα)
    istep = 1 
    io = open_history(fname)
    if !dobo
        Δm = [similar(m) for i in 1:ntries]
    end              
    while true
        if dobo
            idx, foundroot = bostep(m, m0, mnew[istep], χ²[istep], λ²[istep], t, β², λ²GP, F, R, target, lo, hi,
            regularizeupdate=regularizeupdate, ntries=ntries, κ = κ,
            knownvalue=knownvalue, firstvalue=firstvalue, breakonknown=breakonknown)         
        else
            idx, foundroot = occamstep(m, m0, Δm, mnew[istep], χ²[istep], λ²[istep], F, R, target, 
                lo, hi, λ²min, λ²max, β², ntries, knownvalue=knownvalue, regularizeupdate = regularizeupdate)
        end
        prefix = isempty(fname) ? fname : fname*" : "
        @info prefix*"iteration: $istep χ²: $(χ²[istep][idx]) target: $target"
        m = mnew[istep][idx]
        oidx[istep] = idx
        isa(io, Nothing) || write_history(io, [istep; χ²[istep][idx]/target₀; vec(m)])
        foundroot && break
        noimprovement = iszero(λ²[istep][idx][2])  ? true : false
        if (istep == nstepsmax - 1) || noimprovement 
            target = χ²[istep][idx] # exit with smoothest
        end    
        istep += 1
        istep > nstepsmax && break
    end
    isa(io, Nothing) || begin @info "Finished "*fname; close(io) end
    return map(x->x[1:istep], (mnew, χ², λ², oidx))
end    

# nuisance inversion stuff
function gradientinv(   m::AbstractVector,
                        m0::AbstractVector,
                        nu::AbstractVector, 
                        F::Operator; 
                        regtype=:R1,
                        nstepsmax = 10,
                        ntries = 6,
                        target = nothing,
                        lo = -3.,
                        hi = 1.,
                        λ²min = 0,
                        λ²max = 8,
                        β² = 0.,
                        nuλ²frac = zeros(0),
                        nubounds = zeros(0),
                        ndivsnu = zeros(Int, 0),
                        ntriesnu = 10,
                        regularizeupdate = false,
                        knownvalue=0.7,
                        firstvalue=:middle,
                        κ = GP.Mat52(),
                        acqfun = GP.EI(),
                        δtry = 1e-2,
                        breakonknown=true,
                        breaknuonknown=false,
                        reducenuto=0.2,
                        fname="")
    R = makereg(regtype, F)                
    ndata = length(F.res)
    isnothing(target) && (target = ndata)
    target₀ = target
    mnew = [[similar(m) for i in 1:ntries] for j in 1:nstepsmax]
    χ²   = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
    χ²nu   = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
    λ² = [Vector{Vector{Float64}}(undef, 0) for j in 1:nstepsmax]
    nunew = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
    Δm = [similar(m) for i in 1:ntries]
    oidx = zeros(Int, nstepsmax)
    nu_oidx = similar(oidx)  
    ndata = length(F.res)
    t, λ²GP = nuinitbo(nubounds, ndivsnu, nuλ²frac)
    G = GP.KernelStruct(κ, ntriesnu, λ²GP, δtry, t)
    istep = 1 
    io = open_history(fname)
    while true
        # first get optimal nuisance index
        # nu_idx = bostepnu(G, nu, nunew[istep], m, χ²nu[istep], t, λ²GP, F; acqfun, 
        #     ntries=ntriesnu, knownvalue, firstvalue, breakonknown)
        
        # Optim stuff
        f(x) = 2*get_misfit(m, x, F, nubounds)
        f_abstol = breaknuonknown ? reducenuto*f(nu) : 0.  
        nu = Optim.minimizer(optimize(f, nu, BFGS(), 
            Optim.Options(show_trace=true, f_abstol=f_abstol, iterations=ntriesnu, successive_f_tol=0)))
        nunew[istep] = nu
        #
        idx, foundroot = occamstep(m, m0, Δm, mnew[istep], χ²[istep], λ²[istep], F, R, target, 
            lo, hi, λ²min, λ²max, β², ntries, knownvalue=knownvalue, regularizeupdate = regularizeupdate)
        prefix = isempty(fname) ? fname : fname*" : "
        @info prefix*"iteration: $istep χ²: $(χ²[istep][idx]) target: $target"
        m = mnew[istep][idx]
        # nu = nunew[istep][nu_idx]
        oidx[istep] = idx
        # nu_oidx[istep] = nu_idx
        isa(io, Nothing) || write_history(io, [istep; χ²[istep][idx]/target₀; vec(m)])
        foundroot && break
        noimprovement = iszero(λ²[istep][idx][2])  ? true : false
        if (istep == nstepsmax - 1) || noimprovement 
            target = χ²[istep][idx] # exit with smoothest
        end    
        istep += 1
        istep > nstepsmax && break
    end
    isa(io, Nothing) || begin @info "Finished "*fname; close(io) end
    return map(x->x[1:istep], (mnew, nunew, χ², χ²nu, λ², oidx))
end

function nuinitbo(nubounds::Array{Float64, 2}, ndivsnu::Vector{Int}, nufrac::Vector)
    @assert size(nubounds, 1) == length(ndivsnu)
    @assert size(nubounds, 1) == length(nufrac)
    @assert length(nufrac) >= 1. # this is what we divide by to get a fraction of the nuisance range 
    nuranges = map((nu, ndiv)->range(extrema(nu)..., ndiv+1), (eachrow(nubounds)), ndivsnu) # test range for surrogate
    λ²GP =  (abs.(diff(nubounds, dims=2)[:])./nufrac).^2 # length scale square for surrogate
    t = ndgrid(nuranges...)
    tgrid = reduce(vcat, [tt[:]' for tt in t])
    tgrid, λ²GP
end    

function bostepnu(G:: GP.KernelStruct, nu::AbstractVector, nunew::Vector{Vector{Float64}}, m::AbstractVector, χ²::Vector{Float64},
                   t::Array{Float64, 2}, λ²GP::Array{Float64, 1}, F::Operator;

                   ## GP stuff
                   demean = true, 
                   acqfun = GP.EI(),
                   ntries = 10,
                   firstvalue = :middle,
                   highval = 1e32,
                   knownvalue = NaN,
                   breakonknown = false)
    
    χ²₀ = 2*get_misfit(m, nu, F)
    push!(χ², χ²₀) # make sure we get starting value to train
    push!(nunew, nu) # make sure we store starting value
    knownvalue *= χ²₀               
    ttrain = zeros(size(t,1), 0)
    ttrain = hcat(ttrain, nu) # add first training location
    GP.GPfit(G, ttrain, χ²', t)
    # @info "at start of BO: " ttrain, nunew, χ²
    for i = 2:ntries
        nextpos, = getBOsample(G, χ²', ttrain, t, knownvalue, firstvalue, acqfun)
        push!(nunew, t[:, nextpos])
        temp = 2*get_misfit(m, nunew[i], F)
        temp = isnan(temp) ? highval : temp
        push!(χ², temp)                     # next training value
        ttrain = hcat(ttrain, t[:,nextpos]) # next training location
        # @info "at BO try $i" ttrain, χ², nunew
        (χ²[i] <= knownvalue && breakonknown) && break
    end
    idx = argmin(χ²)
    get_misfit(m, nunew[idx], F) # make sure that geometry is updated
    idx
end

function open_history(fname)
    @assert !isfile(fname) "$fname exists"
    if !isempty(fname)
        io = open(fname, "w")
    end
end        

function write_history(io, v::Vector)
    for el in v
        msg = @sprintf("%e\t", el)
        write(io, msg)
    end
    write(io, "\n")
end

# L1 experimental stuff
function coordinatedesc(A, x::AbstractVector, y::AbstractVector, λ², W)
    a = getaⱼ(A)
    for l in λ²
        for j in 1:size(A, 2)
            cⱼ = getcⱼ(A, x, y, j)
            x[j] = getx̂(a[j], cⱼ, l)
        end
        r = (y-A*x)./W
        ϕd = sum(r'r/length(y))
        ϕd <= 1. && (@info l; break)
    end    
end

function getaⱼ(A)
    2*[sum(A[:,j].^2) for j in 1:size(A,2)]
end    

function getcⱼ(A, x::AbstractVector, y::AbstractVector, j::Int)
    cⱼ = 0.
    for i = 1:lastindex(y)
        cⱼ +=  (dot(A[i,:], x) - A[:,j][i]*x[j] - y[i])*A[:,j][i]
    end   
    2*cⱼ 
end

function getx̂(aⱼ, cⱼ, λ²)
    if λ² < cⱼ
        (-cⱼ + λ²)/aⱼ
    elseif cⱼ < -λ²  
        -(cⱼ + λ²)/aⱼ
    else
        0.
    end          
end

struct Int1Dop
    η
end

struct Int1Dtop
    η
end

Int1Dtop(x::Int1Dop) = Int1Dtop(x.η)

function (foo::Int1Dop)(x::Vector)
    vcat(x[1]/foo.η, cumsum(x)[2:end])
end

function (foo::Int1Dtop)(x::Vector)
    y = reverse(cumsum(reverse(x)[1:end-1]))
    vcat(y[1]+x[1]/foo.η, y)
end    

function makeinverseR1(n; η=1)
    Dinv = Int1Dop(η)
    Dtinv = Int1Dtop(Dinv)
    LinearMap(Dinv, Dtinv, n)
end