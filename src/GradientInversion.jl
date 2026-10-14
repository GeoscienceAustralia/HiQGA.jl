using LinearMaps, SparseArrays, PositiveFactorizations, Lasso, FFTW
using Wavelets
import Roots
using .AbstractOperator, Optim
import .AbstractOperator.setnuforinvtype
import .AbstractOperator.modeldim

abstract type L1Transform end

function makeregR0(F::Operator)
    n = modeldim(F)
    LinearMap(I, n)
end

function makeregR1(F::Operator1D)
    n = modeldim(F)
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

# this reflect capability is not usually used
function pushback_reflect(m, lo, hi; reflect=true)
    for i in eachindex(m)
        while (m[i]<lo) || (m[i]>hi)
                (m[i]<lo) && (m[i] = reflect ? 2*lo - m[i] : lo)
                (m[i]>hi) && (m[i] = reflect ? 2*hi - m[i] : hi)
        end
    end
end  

# push model back within bounds
pushback(m, lo, hi) = clamp!(m, lo, hi)

function newtonstep(m::AbstractVector, m0::AbstractVector, F::Operator, λ²::Float64, β²::Float64, R::SparseMatrixCSC)
    JtW = F.J' * F.W
    Wr  = F.W * F.res
    WJm = F.W * (F.J * m)
    H   = (JtW * JtW' + λ²*R'R + λ²*β²*I)
    U   = cholesky(Positive, H, Val{false}).U
    # Solves for the new model m directly. Equivalent to
    # H·m_new = JᵀW·(W·J·m_lin − W·res) + λ²β²·m₀, i.e., regularises the model itself
    # (not its update).
    U \ (U' \ (JtW * (WJm .- Wr) .+ λ² * β² * m0))
end

function cd_lasso_linear_substep(A, C, y, m, m0, l², β²; cd_maxiter, cd_tol,
                               penalty_factor::Vector{Float64})
    ndata = size(A, 1)
    # penalty_factor: per-coefficient L1 weight ωⱼ (penalty becomes λ²ωⱼ on coeff j).
    # ωⱼ=0 frees coeff j (no shrinkage). standardizeω=false keeps the weights
    # literal — otherwise Lasso.jl rescales them to sum to ncoef, shifting λ².
    T = fit(LassoPath, A, y;
            λ = [l²]/ndata,
            standardize = false, intercept = false,
            cd_maxiter, cd_tol,
            maxncoef = size(A, 2),
            penalty_factor = penalty_factor, standardizeω = false)
    z = Matrix(T.coefs)[:]
    return softthresh.(C * z, β² * l²) .+ m0
end

function _occamstep_core(m::AbstractVector, m_proposed::AbstractVector,
                         mnew::Vector{Vector{Float64}},
                         χ²::Vector{Float64}, λ²::Vector{Vector{Float64}},
                         F::Operator, target, lo, hi, λ²min, λ²max, ntries,
                         step_fn, phase2_fn, χ²₀;
                         knownvalue=NaN, noroot=false, countmax=6,
                         sweep_inner_mode::Symbol=:never,
                         sweep_inner_iters::Int=3,
                         sweep_inner_tol::Float64=1e-3,
                         sweep_inner_trigger::Float64=1.5,
                         state_refresh::Union{Function,Nothing}=nothing,
                         sb_reset::Union{Function,Nothing}=nothing,
                         verbose=false)
    knownvalue *= χ²₀
    α = 1
    count = 0
    λ²_seq = 10 .^reverse(LinRange(λ²min, λ²max, ntries))
    while true
        for (i, l²) in enumerate(λ²_seq)
        
            if count == 0
                m_proposed[i] = step_fn(l²)
            end
            # Convex combination of current m and the proposed step (line-search by α).
            @. mnew[i] = α * m_proposed[i] + (1 - α) * m
            pushback(mnew[i], lo, hi)
            chi2 = getχ²(F, mnew[i])
            nu = [l²; α]
            count == 0 && (push!(χ², chi2); push!(λ², nu))
            χ²[i] = chi2
            λ²[i] = nu
            if χ²[i] <= knownvalue || χ²[i] <= target
                if i<ntries
                    map(x->deleteat!(x, i+1:length(x)), (mnew, χ², λ²))
                end
                break
            end
        end
        if all(χ² .>= χ²₀)
            α = count < countmax - 1 ? α/2 : 0.  
        else
            count = countmax
        end
        count += 1
        count > countmax && break
    end
    foundroot = false
    idx = argmin(χ²)
    target_met = any(χ² .< target)
    target_met && (idx = findfirst(χ² .<= target))
    if target_met
        if idx > 1 && !noroot
            a = λ²[idx][1]      # left end of bracketing λ² (reversed in sweep)
            b = λ²[idx-1][1]    # right end on the number line
            phase2_fn(idx, a, b)
        end
        foundroot = true
    end
    idx, foundroot
end

# L2 method: dispatches when caller passes R::SparseMatrixCSC.
# Newton step + Cholesky factorisation.
function occamstep(m::AbstractVector, m0::AbstractVector, m_proposed::AbstractVector,
                   mnew::Vector{Vector{Float64}},
                   χ²::Vector{Float64}, λ²::Vector{Vector{Float64}},
                   F::Operator, R::SparseMatrixCSC, target,
                   lo, hi, λ²min, λ²max, β², ntries;
                   knownvalue=NaN,
                   noroot=false, countmax=6,
                   cd_maxiter=200_000, cd_tol=1e-7, adjtest=false, timeit=false,  # adjtest ignored (Cholesky L2)
                   interproot=false, interptol=1.0, usesecant=true,
                   secant_maxiter=2, force_secant_iters=false,   # accepted but unused (no IRLS in Cholesky L2)
                   phase2_inner_iters=3, phase2_inner_tol=1e-3, phase2_plateau_max_stalls=2,  # accepted but unused
                   sweep_inner_mode::Symbol=:never, sweep_inner_iters::Int=3,           # accepted but unused
                   sweep_inner_tol::Float64=1e-3, sweep_inner_trigger::Float64=1.5,     # accepted but unused
                   verbose=false)                                  # accepted but unused
    χ²₀ = getχ²(F, m, computeJ=true)
    step_fn(l²)          = newtonstep(m, m0, F, l², β², R)
    phase2_fn(idx, a, b) = dophase2_L2(m, m0, F, R, lo, hi,
                                       target, a, b, mnew, χ², λ², β², idx; timeit, interproot, interptol, usesecant)
    _occamstep_core(m, m_proposed, mnew, χ², λ², F, target, lo, hi, λ²min, λ²max, ntries,
                    step_fn, phase2_fn, χ²₀; knownvalue, noroot, countmax)
end

# L1 method: dispatches when caller passes an L1Transform.
# Materialises T = W·J·C once, then drives Lasso.jl at each l².
function occamstep(m::AbstractVector, m0::AbstractVector, m_proposed::AbstractVector,
                   mnew::Vector{Vector{Float64}},
                   χ²::Vector{Float64}, λ²::Vector{Vector{Float64}},
                   F::Operator, t::L1Transform, target,
                   lo, hi, λ²min, λ²max, β², ntries;
                   knownvalue=NaN,
                   noroot=false, countmax=6,
                   cd_maxiter=200_000, cd_tol=1e-7, adjtest=false,  # adjtest ignored (L1)
                   timeit=false, interproot=false, interptol=1.0, usesecant=true,
                   secant_maxiter=2, force_secant_iters=false,   # accepted but unused (no IRLS in L1 path)
                   phase2_inner_iters=3, phase2_inner_tol=1e-3, phase2_plateau_max_stalls=2,  # accepted but unused
                   sweep_inner_mode::Symbol=:never, sweep_inner_iters::Int=3,           # accepted but unused
                   sweep_inner_tol::Float64=1e-3, sweep_inner_trigger::Float64=1.5,     # accepted but unused
                   verbose=false)                                  # accepted but unused
    χ²₀  = getχ²(F, m, computeJ=true)                   
    C    = makeC(t)
    mat_timed = @timed Matrix(F.W*F.J*C)
    A    = mat_timed.value
    timeit && @info "materialise A = WJC: $(round(mat_timed.time, sigdigits=4)) s " *
                    "($(size(A,1))×$(size(A,2)), once per outer iteration)"
    y    = -F.W*(F.res - F.J*(m-m0))
    penalty_factor = makepenaltyfactor(t, size(A, 2))   # per-coefficient L1 weights (weight[1]=0 by convention)
    function step_fn(l²)
        try
            cd_lasso_linear_substep(A, C, y, m, m0, l², β²; cd_maxiter, cd_tol, penalty_factor)
        catch e
            @info "whoops!!! main body linearisation on pid $(myid())"
            @info e
            zeros(size(m))
        end
    end
    
    phase2_fn(idx, a, b) = dophase2_L1(A, C, y, m, m0, F, β², lo, hi,
                                       target, a, b, mnew, χ², λ², idx;
                                       cd_maxiter, cd_tol, penalty_factor, timeit, interproot, interptol, usesecant)
    _occamstep_core(m, m_proposed, mnew, χ², λ², F, target, lo, hi, λ²min, λ²max, ntries,
                    step_fn, phase2_fn, χ²₀; knownvalue, noroot, countmax)
end

function dophase2_L2(m, m0, F, R, lo, hi, target, a, b, mnew,  χ², λ², β², idx; timeit=false, interproot=false, interptol=1.0, usesecant=true)
    α = λ²[idx][2]
    p2 = @timed begin
        # pure Occam root-find (opt-out, interproot=false)
        try
            λ²[idx][1] = 10 .^Roots.find_zero(log10l² -> fλ²_L2(α, m, m0, F, 10 .^log10l², β², R, lo, hi) - target, (log10(a), log10(b)))
            mp = newtonstep(m, m0, F, λ²[idx][1], β², R)
            @. mnew[idx] = α * mp + (1 - α) * m
            pushback(mnew[idx], lo, hi)
            χ²[idx] = getχ²(F, mnew[idx])
        catch e
            @info "L2 root-finding skipped on pid $(myid()) (bracket invalid or numerical issue), retaining sweep model at idx"
            @info e
            # leave mnew[idx]/λ²[idx]/χ²[idx] at sweep values
        end
    end
    timeit && @info "L2 (Cholesky) phase 2 ($(interproot ? "interp" : "root-find")): $(round(p2.time, sigdigits=4)) s"
    nothing
end    

function dophase2_L1(A, C, y, m, m0, F, β², lo, hi, target, a, b, mnew, χ², λ², idx; cd_maxiter=200_000, cd_tol=1e-7, penalty_factor::Vector{Float64}, timeit=false, interproot=false, interptol=1.0, usesecant=true)
    α = λ²[idx][2]
    p2 = @timed begin
        # pure Occam root-find (opt-out, interproot=false)
        try
            λ²_root = 10 .^Roots.find_zero(log10l² -> fλ²_L1(A, C, y, α, m, m0, F, 10 .^log10l², β², lo, hi; cd_maxiter, cd_tol, penalty_factor) - target, (log10(a), log10(b)), Roots.Brent())
            m_proposed = cd_lasso_linear_substep(A, C, y, m, m0, λ²_root, β²; cd_maxiter, cd_tol, penalty_factor)
            # both root-finding and the final substep succeeded — commit
            λ²[idx][1] = λ²_root
            @. mnew[idx] = α * m_proposed + (1 - α) * m
            pushback(mnew[idx], lo, hi)
            χ²[idx] = getχ²(F, mnew[idx])
        catch e
            @info "L1 root-finding skipped on pid $(myid()) (bracket invalid or numerical issue), retaining line-search idx"
            @info e
            # leave mnew[idx]/λ²[idx]/χ²[idx] at line-search values
        end
    end
    timeit && @info "L1 phase 2 ($(interproot ? "interp" : "root-find")): $(round(p2.time, sigdigits=4)) s"
    nothing
end 

function fλ²_L2(α, m, m0, F, l², β², R, lo, hi)
    m_proposed = newtonstep(m, m0, F, l², β², R)
    getnewχ²(m, α, m_proposed, lo, hi, F)
end

function fλ²_L1(A, C, y, α, m, m0, F, l², β², lo, hi; cd_maxiter=200_000, cd_tol=1e-7, penalty_factor::Vector{Float64})
    try
        m_proposed = cd_lasso_linear_substep(A, C, y, m, m0, l², β²; cd_maxiter, cd_tol, penalty_factor)
        getnewχ²(m, α, m_proposed, lo, hi, F)
    catch e
        @info "whoops!!! root finding beginning on pid $(myid())"
        @info e
        getχ²(F, m)
    end
end

function getnewχ²(m, α, m_proposed, lo, hi, F)
    mnew = @. α * m_proposed + (1 - α) * m
    pushback(mnew, lo, hi)                
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

function getχ²(F, m, mn; computeJ=false)
    r, W = F.res, F.W
    r₀ = copy(r)    
    getresidual(F, m, mn)
    χ² = norm(W*r)^2
    computeJ || (r .= r₀)
    return χ²
end

# no nuisances version
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
                        β² = 0.,
                        knownvalue=0.7,
                        verbose = true,
                        minimprovfrac = nothing,
                        minimprovkickinstep = round(Int, nstepsmax/2),
                        noroot = false,
                        countmax = 6,
                        cd_maxiter = 200_000, 
                        cd_tol = 1e-6,
                        L1transform::Union{L1Transform,Nothing} = nothing,
                        l2op = nothing,
                        lassoirlsop = nothing,
                        lassosbop = nothing,
                        adjtest=false,
                        timeit=false,
                        interproot=false,
                        interptol=1.0,
                        usesecant=false,
                        secant_maxiter=2,
                        force_secant_iters=false,
                        phase2_inner_iters=3,
                        phase2_inner_tol=1e-3,
                        phase2_plateau_max_stalls=2,
                        sweep_inner_mode=:never,
                        sweep_inner_iters=3,
                        sweep_inner_tol=1e-3,
                        sweep_inner_trigger=1.5,
                        fname="")
    countmax < 1 && @warn("countmax must be at least 1 to guarantee no uphill")
    # Dispatch priority (most specific first):
    # L1transform → L1-CD route; else sparse R → Cholesky L2.
    reg = (L1transform === nothing ? makereg(regtype, F) : L1transform)
    ndata = length(F.res)
    isnothing(target) && (target = ndata)
    target₀ = target
    mnew = [[similar(m) for i in 1:ntries] for j in 1:nstepsmax]
    χ²   = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
    λ² = [Vector{Vector{Float64}}(undef, 0) for j in 1:nstepsmax]
    oidx = zeros(Int, nstepsmax)  
    ndata = length(F.res)
    istep = 1 
    io = open_history(fname)
    m_proposed = [similar(m) for i in 1:ntries]
    littleimprovement = false              
    while true
        # (reg dispatches: SparseMatrixCSC→L2 Cholesky, L1Transform→L1-CD)
        idx, foundroot = occamstep(m, m0, m_proposed, mnew[istep], χ²[istep], λ²[istep], F, reg, target, 
            lo, hi, λ²min, λ²max, β², ntries; 
            knownvalue, noroot, countmax, cd_maxiter, cd_tol, adjtest, timeit, interproot, interptol, usesecant, secant_maxiter, force_secant_iters, phase2_inner_iters, phase2_inner_tol, phase2_plateau_max_stalls, sweep_inner_mode, sweep_inner_iters, sweep_inner_tol, sweep_inner_trigger, verbose)
        prefix = isempty(fname) ? fname : fname*" : "
        verbose && @info prefix*"iteration: $istep χ²: $(χ²[istep][idx]) target: $target"
        m = mnew[istep][idx]
        oidx[istep] = idx
        isa(io, Nothing) || write_history(io, [istep; χ²[istep][idx]/target₀; vec(m)])
        foundroot && break
        prevχ² = istep > 1 ? χ²[istep-1][oidx[istep-1]] : 1e9
        if !isnothing(minimprovfrac) && (istep > minimprovkickinstep)
            if (prevχ² - χ²[istep][idx])/prevχ² < minimprovfrac
                littleimprovement = true
            end
        end
        noimprovement = iszero(λ²[istep][idx][2])  ? true : false
        if istep > 1
            # check if previous is also stuck
            iszero(λ²[istep-1][oidx[istep-1]][2]) && noimprovement && break
        end
        if (istep == nstepsmax - 1) || noimprovement || littleimprovement
            target = χ²[istep][idx] # exit with smoothest if possible
        end
        # check for coordinate descent fail -- not sure if this is necessary for ordinary L2 
        # α not zero so step size is not zero so noimprovement is false
        !noimprovement && (χ²[istep][idx] == prevχ²) && break    
        istep += 1
        istep > nstepsmax && break
    end
    if !isnothing(io) 
        verbose && @info "Finished "*fname
        close(io) 
        nothing
        return
    end
    # sometimes a smoothest root is not found
    (istep > nstepsmax) && (istep -= 1)
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
                        nubounds = zeros(0),
                        ntriesnu = 10,
                        knownvalue=0.7,
                        breaknuonknown=false,
                        reducenuto=0.2,
                        debuglevel = 0,
                        usebox = true,
                        boxiters = 2,
                        verbose = true,
                        minimprovfrac = nothing,
                        minimprovkickinstep = round(Int, nstepsmax/2),
                        noroot = false,
                        countmax = 6,
                        cd_maxiter = 200_000, 
                        cd_tol = 1e-6,
                        L1transform::Union{L1Transform,Nothing} = nothing,
                        l2op = nothing,
                        lassoirlsop = nothing,
                        lassosbop = nothing,
                        adjtest=false,
                        timeit=false,
                        interproot=false,
                        interptol=1.0,
                        usesecant=false,
                        secant_maxiter=2,
                        force_secant_iters=false,
                        phase2_inner_iters=3,
                        phase2_inner_tol=1e-3,
                        phase2_plateau_max_stalls=2,
                        sweep_inner_mode=:never,
                        sweep_inner_iters=3,
                        sweep_inner_tol=1e-3,
                        sweep_inner_trigger=1.5,
                        fname="")
    # Dispatch priority (most specific first):
    # L1transform → L1-CD route; else sparse R → Cholesky L2.
    reg = (L1transform === nothing ? makereg(regtype, F) : L1transform)
    ndata = length(F.res)
    isnothing(target) && (target = ndata)
    target₀ = target
    mnew = [[similar(m) for i in 1:ntries] for j in 1:nstepsmax]
    χ²   = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
    χ²nu   = Vector{Float64}(undef, nstepsmax)
    λ² = [Vector{Vector{Float64}}(undef, 0) for j in 1:nstepsmax]
    nunew = [Vector{Float64}(undef, 0) for j in 1:nstepsmax]
    m_proposed = [similar(m) for i in 1:ntries]
    oidx = zeros(Int, nstepsmax)
    ndata = length(F.res)
    istep = 1 
    io = open_history(fname)
    littleimprovement = false  
    while true
        # Optim stuff for nuisances
        # f(x) = 2*get_misfit(m, x, F, nubounds) # why? always set usebox=true
        f(x) = getχ²(F, m, setnuforinvtype(F, x)) 
        f_abstol = breaknuonknown ? reducenuto*f(nu) : 0.
        show_trace = debuglevel > 0 ? true : false
        if usebox
            res = optimize(f, nubounds[:,1], nubounds[:,2], nu, Fminbox(BFGS(linesearch=Optim.LineSearches.BackTracking(order=3))), 
                Optim.Options(;show_trace, outer_f_abstol=f_abstol, f_abstol=f_abstol, outer_iterations = boxiters, iterations=ntriesnu)) 
        else
            @warn "will not be stable if nuisance outside bounds! set usebox=true"   
            res = optimize(f, nu, BFGS(), 
                Optim.Options(;show_trace, f_abstol, iterations=ntriesnu, successive_f_tol=0))
        end        
        (debuglevel > 1) && @info res
        nu = res.minimizer    
        nunew[istep] = nu
        χ²nu[istep] = res.minimum
        # normal Occam for conductivities (reg dispatches: SparseMatrixCSC→L2 Cholesky, L1Transform→L1-CD)
        idx, foundroot = occamstep(m, m0, m_proposed, mnew[istep], χ²[istep], λ²[istep], F, reg, target, 
            lo, hi, λ²min, λ²max, β², ntries; knownvalue, noroot, countmax, cd_maxiter, cd_tol, adjtest, timeit, interproot, interptol, usesecant, secant_maxiter, force_secant_iters, phase2_inner_iters, phase2_inner_tol, phase2_plateau_max_stalls, sweep_inner_mode, sweep_inner_iters, sweep_inner_tol, sweep_inner_trigger, verbose)
        prefix = isempty(fname) ? fname : fname*" : "
        verbose && @info prefix*"iteration: $istep χ²: $(χ²[istep][idx]) target: $target"
        m = mnew[istep][idx]
        oidx[istep] = idx
        isa(io, Nothing) || write_history(io, [istep; χ²[istep][idx]/target₀; vec(m); nu])
        foundroot && break
        prevχ² = istep > 1 ? χ²[istep-1][oidx[istep-1]] : 1e9
        if !isnothing(minimprovfrac) && (istep > minimprovkickinstep)
            if (prevχ² - χ²[istep][idx])/prevχ² < minimprovfrac
                littleimprovement = true
            end
        end
        noimprovement = iszero(λ²[istep][idx][2])  ? true : false
        if istep > 1
            # check if previous is also stuck
            iszero(λ²[istep-1][oidx[istep-1]][2]) && noimprovement && break
        end
        if (istep == nstepsmax - 1) || noimprovement || littleimprovement
            target = χ²[istep][idx] # exit with smoothest if possible
        end
        # check for coordinate descent fail -- not sure if this is necessary for ordinary L2 
        # α not zero so step size is not zero so noimprovement is false
        !noimprovement && (χ²[istep][idx] == prevχ²) && break
        istep += 1
        istep > nstepsmax && break
    end
    if !isnothing(io) 
        verbose && @info "Finished "*fname
        close(io) 
        nothing
        return
    end
    # nuisances mess up occam so this is a kludge as 
    # sometimes a smoothest root is not found
    (istep > nstepsmax) && (istep -= 1)
    return map(x->x[1:istep], (mnew, nunew, χ², χ²nu, λ², oidx))
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
    flush(io)
end

# L1 experimental stuff
function coordinatedesc(A, x::AbstractVector, y::AbstractVector, λ², W; tol=1e-3, maxiters=1_000, verbose=false)
    aⱼ = getaⱼ(A)
    lopt = 0
    for l in λ²
        cdstep!(A, aⱼ, x, y, l; tol, maxiters, verbose)
        r = (y-A*x)./W
        ϕd = r'r/length(y)
        ϕd <= 1. && (@info l, r'r; lopt=l; break)
    end
    lopt    
end

function getaⱼ(A)
    2*[sum(A[:,j].^2) for j in 1:size(A,2)]
end    

function cdstep!(A, aⱼ, x, y, l; tol=1e-3, maxiters=1_000, verbose=false)
    maxdel = 0.
    maxdel_prev = 0.
    converged = false
    for iter in 1:maxiters
        maxdel_prev = maxdel
        for j in 1:size(A, 2)
            cⱼ = getcⱼ(A, x, y, j)
            xjold = x[j]
            x[j] = getx̂(aⱼ[j], cⱼ, l)
            maxdel = max(maxdel, abs.(xjold-x[j]))
        end
        if iter>=2
            diffdelfrac = abs((maxdel-maxdel_prev)/maxdel)
            verbose && @info "iter: $iter out of $maxiters, maxdelta: $maxdel, change in delta: $diffdelfrac"
            if (diffdelfrac < tol) | isnan(diffdelfrac) # less than tol, or maxdel is 0.
                converged = true
                break
            end
        end 
    end
    if !converged 
        @error "did not converge"
    end 
    return nothing   
end    

function getcⱼ(A, x::AbstractVector, y::AbstractVector, j::Int)
    cⱼ = 0.
    for i = 1:lastindex(y)
        cⱼ +=  (dot(A[i,:], x) - A[:,j][i]*x[j] - y[i])*A[:,j][i]
    end   
    2*cⱼ 
end

function getx̂(aⱼ, cⱼ, λ²)
    # if λ² < cⱼ
    #     β = (-cⱼ + λ²)/aⱼ
    # elseif cⱼ < -λ²  
    #     β = -(cⱼ + λ²)/aⱼ
    # else
    #     β = 0.
    # end
    # which is the same as
    softthresh(-cⱼ/aⱼ, λ²/aⱼ)
    # β
end

function softthresh(β, λ)
    # sign(β)*max(0, abs(β)-λ)
    # more stably
    if β < -λ
        β = β + λ
    elseif β > λ  
        β = β - λ
    else
        β = 0.
    end
    β 
end

# The logic should be for Int1Dop{T}
# freelevel true -> weights are 0 on w[1] (zero penalty on reg side)
# free_level false -> z[1]/η is 0, doesn't matter what w[1] is (operator side weighting)

struct Int1Dop{T}
    η::T
    free_level::Bool
end
Int1Dop(η) = Int1Dop(η, true)   # backward-compat default: keep the η^{-1} level column

struct Int1Dtop{T}
    η::T
    free_level::Bool
end
Int1Dtop(η) = Int1Dtop(η, true)

Int1Dtop(x::Int1Dop) = Int1Dtop(x.η, x.free_level)

# Forward action. When free_level=false the first column of the operator matrix
# is identically zero — the level-carrying η^{-1} entry is zeroed, so z[1] does
# not contribute to the output. Same operator shape (p × p), same code path.
function (foo::Int1Dop)(x::AbstractVector)
    t = foo.free_level ? x[1]/foo.η : zero(eltype(x))
    cumsum([t, x[2:end]...])
end

# Adjoint action. When free_level=false the first row of the adjoint (i.e., the
# transposed level column) is identically zero.
function (foo::Int1Dtop)(x::AbstractVector)
    t = cumsum(reverse(x))
    if foo.free_level
        t[end] /= foo.η
    else
        t[end] = zero(eltype(x))
    end
    reverse(t)
end

function makeinverseR1(n; η=1, free_level=true)
    Dinv  = Int1Dop(η, free_level)
    Dtinv = Int1Dtop(Dinv)
    LinearMap(Dinv, Dtinv, n)
end

struct R1Transform <: L1Transform
    n::Int
    η::Float64
    free_level::Bool   # if true, coefficient 1 (the model level) is unpenalised (Λ₁=0)
    weights::Vector{Float64}   # per-coefficient L1 penalty weights; ω[1]=0 iff free_level
    function R1Transform(n::Int, η::Real, free_level::Bool)
        w = ones(n)
        free_level && (w[1] = 0.0)
        new(n, Float64(η), free_level, w)
    end
end
# η=1 ⟹ makeinverseR1 is exact cumsum (well-conditioned, bijective). The level is
# freed on the PENALTY side via free_level (Λ₁=0), not by shrinking η — see §5.1 of
# the manuscript.
R1Transform(F::Operator; η=1, free_level=true) = R1Transform(modeldim(F), η, free_level)
R1Transform(n::Int; η=1, free_level=true) = R1Transform(n, η, free_level)

struct R2Transform <: L1Transform
    n::Int
    η::Float64
    weights::Vector{Float64}   # ω[1]=ω[2]=0 (double integrator has two constants of integration)
    function R2Transform(n::Int, η::Real)
        w = ones(n); w[1] = 0.0; w[2] = 0.0
        new(n, Float64(η), w)
    end
end
R2Transform(F::Operator; η=1e-3) = R2Transform(modeldim(F), η)
R2Transform(n::Int; η=1e-3) = R2Transform(n, η)

# Separable N-D transform. `factors` are stored in math/kron order, so the
# transform is kron(C_factors[1], C_factors[2], ...) = C_factors[1] ⊗ C_factors[2] ⊗ ...
# i.e. factor[1] is the leftmost (slowest-varying) Kronecker factor.
#
# Two evaluation paths, selected by `fwdonly`:
#  - fwdonly=false (default): adjoint path via LinearMaps' kron(...). Handles ANY
#    number of factors, but each factor's makeC must define a valid adjoint.
#    Use for differences (R1/R2) and ORTHOGONAL wavelets.
#  - fwdonly=true: adjoint-free path via kronmul_fwdonly_vec. Currently 2-factor
#    only. REQUIRED for biorthogonal wavelets (CDF etc.) where dwt ≠ idwtᵀ.
#
# Footgun note: the adjoint path with a biorthogonal wavelet factor silently gives
# wrong results (the dwt adjoint is invalid). For biorthogonal wavelets always set
# fwdonly=true. For orthogonal wavelets / differences either path works and must
# agree bit-for-bit — a built-in cross-check.
struct KronTransform <: L1Transform
    factors::Vector{<:L1Transform}
    fwdonly::Bool
end
KronTransform(factors::Vector{<:L1Transform}; fwdonly=false) = KronTransform(factors, fwdonly)

# Identity factor: a dispatch shim so LinearAlgebra.I can sit in a KronTransform's
# typed factor list. makeC returns LinearMap(I, n); it carries no penalty. Lets the
# directional atoms be written as their Kronecker formulae, e.g. G_x = Σ_x ⊗ I.
struct IdentityTransform <: L1Transform
    n::Int
    weights::Vector{Float64}   # always ones; no coarse-level DoF of its own
    IdentityTransform(n::Int) = new(n, ones(n))
end

# Directional integrator atoms for 2D anisotropic TV, written as their Kronecker
# formulae (column-major vec, p1=rows=dir 1=y, p2=cols=dir 2=x):
#   G_y = I_{p2} ⊗ Σ_y(p1)  — integrates down rows (direction 1)
#   G_x = Σ_x(p2) ⊗ I_{p1}  — integrates across cols (direction 2); on vec(Z) this is
#                             vec(Z·Σ_xᵀ), realised by LinearMaps' kron. η=1 ⇒ Σ=cumsum.
# Both use fwdonly=false: Σ (the cumsum integrator) and I have valid adjoints, so the
# N-factor adjoint-kron path applies, and materialisation uses only the forward action.
Gy_atom(p1, p2; η=1) = KronTransform([IdentityTransform(p2), R1Transform(p1; η, free_level=false)])
Gx_atom(p1, p2; η=1) = KronTransform([R1Transform(p2; η, free_level=false), IdentityTransform(p1)])

struct DictionaryTransform <: L1Transform
    atoms::Vector{<:L1Transform}
    free_const::Bool   # append an all-ones column, unpenalised (single shared level)
end
DictionaryTransform(atoms::Vector{<:L1Transform}; free_const=false) =
    DictionaryTransform(atoms, free_const)

struct WaveletTransform <: L1Transform
    n::Int
    wt
    weights::Vector{Float64}    # per-coefficient L1 penalty weights, length n
    function WaveletTransform(n::Int, wt, weights::Vector{Float64})
        length(weights) == n ||
            error("WaveletTransform weights length $(length(weights)) ≠ n $n")
        ispow2(n) ||
            error("WaveletTransform n=$n must be a power of 2 (full-depth DWT reduces DC to index 1 only for 2^k lengths)")
        if weights[1] != 0
            @warn "WaveletTransform: weight[1] = $(weights[1]) is nonzero. The DC (scaling) coefficient is unpenalised by convention in Occam-style inversion. If this is intentional (e.g., sparse recovery with known mean), ignore this warning."
        end
        new(n, wt, weights)
    end
end
# default constructors: weight[1]=0 by convention (scaling/DC coefficient is unpenalised),
# other weights default to 1 (unweighted lasso on the detail/high-frequency coefficients).
WaveletTransform(n::Int, wt)              = WaveletTransform(n, wt, _default_l1_weights(n))
WaveletTransform(F::Operator, wt)         = WaveletTransform(modeldim(F), wt)
WaveletTransform(F::Operator, wt, w)      = WaveletTransform(modeldim(F), wt, w)

# Helper: default L1 penalty weights vector — ones with weight[1]=0. Used by the
# convenience constructors of WaveletTransform, CosineTransform, CosineTransform2D
# and SquareWaveletTransform2D so their defaults do not trigger the inner-
# constructor "weight[1] overridden" warning.
function _default_l1_weights(n::Int)
    w = ones(n)
    w[1] = 0.0
    w
end

# Per-coefficient L1 penalty weights (Lasso.jl penalty_factor for CD path; static
# prior Λ_static for IRLS/SB paths). Returns a concrete Vector{Float64} of length
# ncoef. By convention every atom has weight[1] = 0 for its coarse-level DoF
# (model level / DC / scaling / shared constant), OR is used as a factor inside
# a composite (KronTransform / DictionaryTransform) where the freeing lives
# elsewhere — in which case its own weights are all ones.
#   ωⱼ = 1  : penalised at the shared λ²
#   ωⱼ > 1  : suppressed more strongly
#   ωⱼ = 0  : free (no shrinkage) — used for the model level, DC coefficient,
#             wavelet scaling coefficient, and DictionaryTransform's free constant.
#
# GENERIC ACCESSOR — for non-composite transforms whose weights are set at
# construction time and stored in `weights::Vector{Float64}`. Covers R1Transform,
# R2Transform, IdentityTransform, WaveletTransform, SquareWaveletTransform2D,
# CosineTransform, CosineTransform2D. The construction-time invariant (weight[1]=0
# where a free level is present, or ones(n) where none is) is trusted here.
makepenaltyfactor(t::L1Transform, ncoef::Int) = t.weights

# COMPOSITE — KronTransform: the stacked-coefficient penalty is the Kronecker
# product of the factor penalty vectors, in the SAME coefficient order as
# C = kron(makeC.(factors)...). Each factor returns a concrete vector; the per-
# factor weight for 2D coefficient (i,j) is wx[i]·wy[j] — separable Λ = Λ₂⊗Λ₁.
function makepenaltyfactor(t::KronTransform, ncoef::Int)
    factor_ncoef = [size(makeC(f), 2) for f in t.factors]
    ωs = [makepenaltyfactor(f, nc) for (f, nc) in zip(t.factors, factor_ncoef)]
    ω = reduce(kron, ωs)
    length(ω) == ncoef ||
        error("KronTransform penalty length $(length(ω)) ≠ ncoef $ncoef")
    ω
end

# COMPOSITE — DictionaryTransform: concatenate each atom's per-block penalty; if
# there is a free constant column, append one 0 for it (unpenalised, single scalar
# coefficient).
function makepenaltyfactor(t::DictionaryTransform, ncoef::Int)
    atom_ncoef = [size(makeC(a), 2) for a in t.atoms]
    ωs = [makepenaltyfactor(a, nc) for (a, nc) in zip(t.atoms, atom_ncoef)]
    ω = reduce(vcat, ωs)
    t.free_const && (ω = vcat(ω, 0.0))   # the all-ones level column is unpenalised
    length(ω) == ncoef ||
        error("DictionaryTransform penalty length $(length(ω)) ≠ ncoef $ncoef")
    ω
end

makeC(t::R1Transform)   = makeinverseR1(t.n; η = t.η, free_level = t.free_level)

makeC(t::IdentityTransform) = LinearMap(I, t.n)

# R2Transform uses invertible cumsum for both integrators; the "two free levels"
# semantics live entirely on the penalty side (weights[1] = weights[2] = 0).
makeC(t::R2Transform)   = makeinverseR1(t.n; η = t.η) * makeinverseR1(t.n; η = t.η)

# Adjoint path: LinearMaps' kron, any number of factors, needs factor adjoints.
# Forward-only path: adjoint-free, 2-factor only (see kronmul_fwdonly_vec below).
function makeC(t::KronTransform)
    Cs = makeC.(t.factors)
    t.fwdonly || return kron(Cs...)                  # adjoint path, N-factor
    # forward-only path: kron(C₁,C₂)·vec(Z) = vec(C₂·Z·C₁ᵀ) = kronmul_fwdonly_vec(C₂, Z, C₁)
    # so B = last factor, A = first factor; Z is (dim C₂ × dim C₁).
    length(Cs) == 2 || error("fwdonly KronTransform supports exactly 2 factors (got $(length(Cs))); "*
                             "n-factor forward-only not yet implemented — use fwdonly=false for adjoint-having factors")
    A, B = Cs[1], Cs[2]
    nA, nB = size(A, 2), size(B, 2)
    fwd = v -> kronmul_fwdonly_vec(B, reshape(v, nB, nA), A)
    LinearMap(fwd, nA * nB)
end

# Stack atom synthesis maps side by side: C = [C_1 | ... | C_K | (1)]. Each C_i is a
# LinearMap (n × p_i); hcat gives the n × P block-row dictionary. The optional
# all-ones column is the shared model level (freed on the penalty side, see
# makepenaltyfactor). All atoms must map into the same n-dim model space.
#
# ADJOINT SAFETY: occamstep applies this C forward only — Matrix(W*J*C) materialises
# it column-by-column (forward), and C*z synthesises the solution (forward). The
# adjoint is never taken, so a biorthogonal atom (e.g. CDF, whose WaveletTransform
# adjoint dwt≠idwtᵀ) is SAFE to stack here. If a matrix-free path that needs Cᵀ is
# ever added, biorthogonal atoms must be revisited.
function makeC(t::DictionaryTransform)
    Cs = makeC.(t.atoms) # Distpatches to the appropriate method
    n = size(first(Cs), 1) # All atoms must have the same number of rows
    all(C -> size(C, 1) == n, Cs) ||
        error("DictionaryTransform atoms map to differing model dims $(size.(Cs, 1))")
    blocks = t.free_const ? (Cs..., LinearMap(ones(n, 1))) : (Cs...,)
    hcat(blocks...)
end

function makeC(t::WaveletTransform)
    # idwt is the forward synthesis (z → x'); C = idwt is the pure (unweighted)
    # inverse wavelet transform. For ORTHOGONAL wavelets idwt's adjoint is dwt
    # (W⁻¹ = Wᵀ, so idwtᵀ = dwt). Defining the adjoint is what makes
    # kron(makeC(WaveletTransform), ...) work in the 2D adjoint path; for 1D
    # forward-only inversion the adjoint is unused. For a NON-orthogonal (e.g.
    # biorthogonal) wavelet this adjoint is wrong — route through a KronTransform
    # with fwdonly=true (via the WaveletTransform2D constructor), which never
    # calls this adjoint.
    #
    # The scale-dependent taper is NOT baked into the operator any more: it is
    # carried on the penalty side as penalty_factor = t.weights (see
    # makepenaltyfactor below), so C stays the natural orthonormal transform with
    # uniform column scaling, and weights of 0 (free scaling coefficients) are
    # expressible — which the old operator-side C̃ = C·W⁻¹ form could not do.
    fwd = z̃ -> idwt(z̃, t.wt)
    adj = y  -> dwt(y, t.wt)
    LinearMap(fwd, adj, t.n)
end

# Adjoint-free kron action: computes vec(Y) where Y = B·X·Aᵀ (matrix form of
# (A⊗B)·vec(X)) using only forward actions of A and B, via the identity
# M·Aᵀ = (A·Mᵀ)ᵀ. This lets us use kron-like operators even when the constituent
# maps lack a defined adjoint — needed for biorthogonal wavelets (CDF, etc.)
# where dwt is NOT the adjoint of idwt.
#
# Returns vec(Y) directly. The output is column-major vec, which equals reading
# A_BXt row-by-row — implemented via a lazy PermutedDimsArray view to avoid
# materialising the final Y matrix.
function kronmul_fwdonly_vec(B::LinearMap, X::AbstractMatrix, A::LinearMap)
    BX = similar(X, size(B, 1), size(X, 2))
    mul!(BX, B, X)
    BXt = collect(transpose(BX))
    A_BXt = similar(BXt, size(A, 1), size(BXt, 2))
    mul!(A_BXt, A, BXt)
    # vec(Y) where Y = transpose(A_BXt). vec walks the lazy permuted view in
    # column-major order, which equals A_BXt traversed row-by-row → vec(Yᵀᵀ).
    return vec(collect(transpose(A_BXt)))
end

# WaveletTransform2D is now a thin convenience constructor over KronTransform.
# It builds two 1D WaveletTransform factors (x-factor first = math/kron order),
# each carrying its own frequency-taper weights, and routes through the
# forward-only path by default
#
# Anisotropy: pass γy, γx to taper the two directions independently (γ sets both).
# weights default to ones (γ=0 ⇒ unweighted). Factor order [x, y] gives
# kron(Cx, Cy) = Cx ⊗ Cy, matching the validated convention vec(Cy·Z·Cxᵀ).
function WaveletTransform2D(ny::Int, nx::Int, wt; γ=0.0, γy=γ, γx=γ, 
        fwdonly=true, rectangular=true)
    !fwdonly && @warn("ensure that idwt is the adjoint of dwt for $wt")
    if rectangular
        wy = freq_penalty_weights(ny; γ=γy)
        wx = freq_penalty_weights(nx; γ=γx)
        KronTransform([WaveletTransform(nx, wt, wx),    # factor 1: x / columns (math-order leftmost)
                    WaveletTransform(ny, wt, wy)];   # factor 2: y / rows
                    fwdonly)
    else
        
    end
end

# This seems to be the math convention, but inconvenient 
# for level corresponding to wavenumber
function dyadic_level_conventional(i::Int, n::Int, levelmax::Int)
    i <= 2 && return levelmax
    return levelmax - Int(ceil(log2(i / 2)))
end

# This seems better for level corresponding to wavenumber or frequency
function dyadic_level_1d(i::Int, n::Int)
    @assert ispow2(n) "n must be a power of 2 (got $n)"
    levelmax = Int(log2(n))
    @assert i<= 2^levelmax
    -dyadic_level_conventional(i, n, levelmax) + levelmax + 1
end

# higher frequency gets higher penalty weight
function freq_penalty_weights(n::Int; γ=1.0)
    @assert ispow2(n) "n must be a power of 2 (got $n)"
    w = ones(n)
    for i in 1:n
        f = dyadic_level_1d(i, n) - 1
        w[i] = 2.0^(γ * f)
    end
    w[1] = 0.0 # scaling coefficient always zero
    return w
end

# cosine transform stuff AND square wavelet2D
struct CosineTransform <: L1Transform
    n::Int
    weights::Vector{Float64}    # per-coefficient L1 penalty weights, length n
    function CosineTransform(n::Int, weights::Vector{Float64})
        length(weights) == n ||
            error("CosineTransform weights length $(length(weights)) ≠ n $n")
        if weights[1] != 0
            @warn "CosineTransform: weight[1] = $(weights[1]) is nonzero. The DC coefficient is unpenalised by convention in Occam-style inversion. If this is intentional (e.g., sparse recovery with known mean), ignore this warning."
        end
        new(n, weights)
    end
end

# default constructors: weight[1]=0 by convention (DC coefficient is unpenalised),
# other weights default to 1 (unweighted lasso on the detail/high-frequency coefficients).
CosineTransform(n::Int)       = CosineTransform(n, _default_l1_weights(n))
CosineTransform(F::Operator)  = CosineTransform(modeldim(F))
CosineTransform(F::Operator, 
    weights::Vector{Float64}) = CosineTransform(modeldim(F), weights)

struct CosineTransform2D <: L1Transform
    ny::Int
    nx::Int
    weights::Vector{Float64}    # per-coefficient L1 penalty weights, length ny*nx
    function CosineTransform2D(ny::Int, nx::Int, weights::Vector{Float64})
        length(weights) == ny*nx ||
            error("CosineTransform2D weights length $(length(weights)) ≠ ny*nx $(ny*nx)")
        if weights[1] != 0
            @warn "CosineTransform2D: weight[1] = $(weights[1]) is nonzero. The DC coefficient is unpenalised by convention in Occam-style inversion. If this is intentional (e.g., sparse recovery with known mean), ignore this warning."
        end
        new(ny, nx, weights)
    end
end

CosineTransform2D(ny::Int, nx::Int)  = CosineTransform2D(ny, nx, _default_l1_weights(ny*nx))
CosineTransform2D(F::Operator)       = CosineTransform2D(F.ny, F.nx)
CosineTransform2D(F::Operator, weights::Vector{Float64}) =
    CosineTransform2D(F.ny, F.nx, weights)

function makeC(t::CosineTransform)
    # idct is the forward synthesis (z → x'); C = idct is the pure (unweighted)
    # inverse cosine transform. 
    fwd = z̃ -> idct(z̃)
    adj = y  -> dct(y)
    LinearMap(fwd, adj, t.n)
end

function makeC(t::CosineTransform2D)
    # idct is the forward synthesis (z̃ → x'); C = idct is the pure (unweighted)
    # inverse 2-D cosine transform (separable, both axes). collect ⇒ contiguous
    # array so FFTW's stride assertion passes even when handed a view.
    fwd = z̃ -> vec(idct(reshape(collect(z̃), t.ny, t.nx)))
    adj = y  -> vec( dct(reshape(collect(y),  t.ny, t.nx)))
    LinearMap(fwd, adj, t.ny*t.nx)
end

struct SquareWaveletTransform2D <: L1Transform
    n :: Int  # this is already squared, or the total number
    wt
    weights :: Vector{Float64}    # per-coefficient L1 penalty weights, length n
    function SquareWaveletTransform2D(n::Int, wt, weights::Vector{Float64})
        length(weights) == n ||
            error("SquareWaveletTransform2D weights length $(length(weights)) ≠ n $n")
        side = isqrt(n)
        side*side == n ||
            error("SquareWaveletTransform2D n=$n is not a perfect square")
        ispow2(side) ||
            error("SquareWaveletTransform2D side=$side must be a power of 2 (full-depth DWT reduces DC to index 1 only for 2^k sides)")
        if weights[1] != 0
            @warn "SquareWaveletTransform2D: weight[1] = $(weights[1]) is nonzero. The DC (scaling) coefficient is unpenalised by convention in Occam-style inversion. If this is intentional (e.g., sparse recovery with known mean), ignore this warning."
        end
        new(n, wt, weights)
    end
end

SquareWaveletTransform2D(F::Operator, wt) = 
    SquareWaveletTransform2D(modeldim(F), wt, _default_l1_weights(modeldim(F)))
SquareWaveletTransform2D(F::Operator, wt, weights::Vector{Float64}) =
    SquareWaveletTransform2D(modeldim(F), wt, weights)

function makeC(t::SquareWaveletTransform2D)
    side = isqrt(t.n)
    side*side == t.n ||
        error("SquareWaveletTransform2D n=$(t.n) is not a perfect square")
    # idwt on a (side×side) matrix is Wavelets.jl's native 2-D synthesis — the fast
    # path, no Kron. z̃ → reshape → 2-D idwt → vec.
    fwd = z̃ -> vec(idwt(reshape(z̃, side, side), t.wt))
    adj = y  -> vec( dwt(reshape(y,  side, side), t.wt))
    LinearMap(fwd, adj, t.n)
end