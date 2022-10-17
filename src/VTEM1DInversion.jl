module VTEM1DInversion
using ..AbstractOperator, ..CommonToAll
import ..AbstractOperator.get_misfit
import ..main # for McMC
import ..gradientinv # for gradientbased
import ..AbstractOperator.getresidual # for gradientbased
import ..Model, ..Options
using ..AEM_VMD_HMD
import ..AbstractOperator.Sounding # for storing real data
using Random, PyPlot, DelimitedFiles, LinearMaps, SparseArrays, ..GP, LinearAlgebra, Statistics

μ = AEM_VMD_HMD.μ
const pVinv = 1e12

mutable struct dBzdt<:Operator1D
    d          :: Array{Float64, 1}
    useML      :: Bool
    σ          :: Array{Float64, 1}
    F          :: AEM_VMD_HMD.HField
    z          :: Array{Float64, 1}
    nfixed     :: Int
    ρ          :: Array{Float64, 1}
    select     :: Array{Bool, 1}
    ndata      :: Int
    J          :: AbstractArray
    W          :: SparseMatrixCSC
    res        :: Vector
end

function dBzdt(;times           = [1.],
        ramp            = [0. 1],
        rTx             = 10.,  
        zTx             = -30.,
        useML           = false,
        z               = [-1.],
        ρ               = [-1.],
        nfixed          = 1,
        nmax            = length(ρ)+1,
        calcjacobian    = false,
        nfreqsperdecade = 6,
        ntimesperdecade = 10,
        d               = zeros(0),
        σ               = zeros(0),
        showgates       = false
        )
   
    @assert size(σ)  == size(d)
    ndata, select = getndata(d)

    F = AEM_VMD_HMD.HFieldDHT(;
        times,
        ramp,
        nmax,
        zTx,
        rTx,
        rRx         = 0.,
        zRx         = zTx-0.01,
        calcjacobian,
        nfreqsperdecade,
        ntimesperdecade,
    )
    @assert length(F.thickness) >= length(z)
    # for Gauss-Newton
    res, J, W = allocateJ(F.dBzdt_J, σ, select, nfixed, length(ρ), calcjacobian)
    aem = dBzdt(d, useML, σ, F, z, nfixed, copy(ρ), select, ndata, J, W, res)
    showgates && plotwaveformgates(aem)
    aem
end

function allocateJ(FJt, σ, select, nfixed, nmodel, calcjacobian)
    if calcjacobian && !isempty(select)
        J = FJt'
        J = J[select,nfixed+1:nmodel]
        Wdiag = 1 ./σ[select]
        res = similar(Wdiag)
    else    
        res, J, Wdiag = zeros(0), zeros(0), zeros(0)
    end
    W = sparse(diagm(Wdiag))        
    return res, J, W
end

function getndata(d)
    select = .!isnan.(d)
    ndata  = sum(select)
    ndata, select
end

function getresidual(aem::dBzdt, log10σ::Vector{Float64}; computeJ=false)
    F = aem.F
    F.calcjacobian = computeJ
    getfield!(-log10σ, aem)
    f = F.dBzdt[aem.select]
    d = aem.d[aem.select]
    aem.res[:] = f-d
    if computeJ
        select, nfixed = aem.select, aem.nfixed
        copy!(aem.J, F.dBzdt_J'[select, nfixed+1:nfixed+length(log10σ)])
    end
    nothing    
end    

mutable struct VTEMsoundingData <: Sounding
    sounding_string :: String
    X               :: Float64
    Y               :: Float64
    Z               :: Float64
    fid             :: Float64
    linenum         :: Int
    zTx             :: Float64
    rTx             :: Float64
    lowpassfcs      :: Array{Float64, 1}
    times           :: Array{Float64, 1}
    ramp            :: Array{Float64, 2}
    noise           :: Array{Float64, 1}
    data            :: Array{Float64, 1}
end

function VTEMsoundingData(;rRx=nothing, zRx=nothing, zTx=12.,
                            rTx=-12.,lowpassfcs=[],
                            times=[1., 2.], ramp=[1 2; 3 4],
                            noise=[1.], data=[1.], 
                            sounding_string="sounding", X=nothing, Y=nothing, Z=nothing,
                            linenum=nothing, fid=nothing)
    @assert rTx > 0
    @assert zTx < 0 # my coordinate system z down 
    isnothing(rRx) && (rRx = 0.)
    isnothing(zRx) && (zRx = zTx-0.1) # place receiver just above tx centre
    !isempty(lowpassfcs) && @assert all(lowpassfcs .> 0)
    @assert all(diff(times) .>0 )
    @assert all(diff(ramp[:,1]) .>0 )
    @assert all((noise .>0) .| isnan.(noise))
    @assert length(data) == length(noise)
    VTEMsoundingData(sounding_string, X, Y, Z, fid, linenum, zTx, rTx,
        lowpassfcs, times, ramp, noise, data)
end

function read_survey_files(;
    fname_dat="",
    fname_specs_halt="",
    frame_height = -2,
    d        = [-2, -2],
    units=1e-12,
    figsize = (8,4),
    fontsize = 10,
    makesounding = false,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.03,
    X = -1,
    Y = -1,
    Z = -1,
    fid = -1,
    linenum = -1,
    nanchar = "*")

    @assert frame_height > 0
    @assert all(d .> 0)

    @assert X > 0
    @assert Y > 0
    @assert Z > 0
    @assert linenum > 0
    @assert fid > 0
    @info "reading $fname_dat"
    if !isnothing(dotillsounding)
        soundings = readdlm(fname_dat)[startfrom:skipevery:dotillsounding,:]
    else
        soundings = readdlm(fname_dat)[startfrom:skipevery:end,:]
    end
    soundings[soundings .== nanchar] .= NaN
    easting = soundings[:,X]
    northing = soundings[:,Y]
    topo = soundings[:,Z]
    fiducial = soundings[:,fid]
    whichline = soundings[:,linenum]
    d = soundings[:,d[1]:d[2]]
    zTx = -soundings[:,frame_height] # my coordinate system
 
    @info "reading $fname_specs_halt"
    include(fname_specs_halt)
    @assert size(d, 2) == length(times)
    @assert size(d, 2) == length(σ_halt)
    σ_halt[:] .*= units
    d[:]      .*= units
    σ           = sqrt.((multnoise*d).^2 .+ (σ_halt').^2)

    plotsoundingdata(d, σ, times, zTx; figsize, fontsize)
    nsoundings = size(d, 1)
    if makesounding
        s_array = Array{VTEMsoundingData, 1}(undef, nsoundings)
        for is in 1:nsoundings
            l, f = Int(whichline[is]), fiducial[is]
            @info "read $is out of $nsoundings"
            s_array[is] = VTEMsoundingData(;zTx=zTx[is], rTx, 
                times, ramp, noise=σ[is,:], data=d[is,:], 
                sounding_string="sounding_$(l)_$f",
                X=easting[is], Y=northing[is], Z=topo[is], fid=f,
                linenum=l)
        end
        return s_array
    end
end

function plotsoundingdata(d, σ, times, zTx; figsize=(8,4), fontsize=1)
    f, ax = plt.subplots(2, 2, figsize=figsize, gridspec_kw=Dict("width_ratios" => [1,0.01]))
    nsoundings = size(d, 1)
    plot_d = permutedims(d)
    plot_d[plot_d .<0] .= NaN
    img = ax[1].pcolormesh(1:nsoundings, times, log10.(plot_d), shading="nearest")
    ax[1].invert_yaxis()
    cb = colorbar(img, cax=ax[1,2])
    cb.set_label("log10(dBz/dt)")
    ax[1].set_ylabel("time s")
    ax[1].tick_params(labelbottom=false)
    axx = ax[1].twiny()
    axx.semilogy(mean(σ./abs.(d), dims=1)[:], times)
    axx.set_xlabel("avg high alt noise fraction")
    ax[2].plot(1:nsoundings, zTx, label="zTx")
    ax[2].set_xlabel("sounding #")
    ax[2].set_ylabel("height m")
    ax[2].sharex(ax[1])
    ax[2].legend()
    ax[2,2].axis("off")
    nicenup(f, fsize=fontsize)
    nothing
end

# all calling functions here for misfit, field, etc. assume model is in log10 resistivity
# SANS the top. For lower level field calculation use AEM_VMD_HMD structs

function getfield!(m::Model, aem::dBzdt)
    getfield!(m.fstar, aem)
    nothing
end

function getfield!(m::Array{Float64}, aem::dBzdt)
    copyto!(aem.ρ, aem.nfixed+1:length(aem.ρ), 10 .^m, 1:length(m))
    if (aem.ndata>0) | isempty(aem.d)
        AEM_VMD_HMD.getfieldTD!(aem.F,  aem.z, aem.ρ)
    end    
    nothing
end

function get_misfit(m::Model, opt::Options, aem::dBzdt)
    get_misfit(m.fstar, opt, aem)
end

function get_misfit(m::AbstractArray, opt::Options, aem::dBzdt)
    calcmisfit(m, opt.debug, aem)
end

function calcmisfit(m, debug, aem)
    chi2by2 = 0.0
    if !debug
        getfield!(m, aem)
        if aem.ndata>0 
            chi2by2 += getchi2by2(aem.F.dBzdt, aem.d,
                    aem.σ, aem.select, aem.useML, aem.ndata)
        end            
    end
    return chi2by2    
end

function getchi2by2(dBzdt, d, σ, select, useML, ndata)
    r, s, idx = dBzdt, σ, select
    r .= (r - d)./s
    if useML
        chi2by2 = 0.5*ndata*log(norm(r[idx])^2)
    else
        chi2by2 = 0.5*norm(r[idx])^2
    end
end

function computeMLfactor(dBzdt, d, σ, select, ndata)
    if ndata > 0
        r, s, idx = dBzdt, σ, select
        r .= (r - d)./s
        r[idx]'*r[idx]/ndata
    else
        NaN
    end    
end

function computeMLfactor(aem)
    mlfact = computeMLfactor(aem.F.dBzdt, aem.d, aem.σ, aem.select, aem.ndata)
    sqrt(mlfact)
end

# all plotting codes here assume that the model is in log10 resistivity, SANS
# the top layer resistivity. For lower level plotting use AEM_VMD_HMD structs

function plotsoundingcurve(ax, f, t; color=nothing, alpha=1, lw=1)
    if isnothing(color)
        ax.loglog(t, μ*f*pVinv, alpha=alpha, markersize=2, linewidth=lw)
    else
        ax.loglog(t, μ*f*pVinv, color=color, alpha=alpha, markersize=2, linewidth=lw)
    end    
end

function plotdata(ax, d, σ, t; onesigma=true)
    sigma = onesigma ? 1 : 2
    ax.errorbar(t, μ*d*pVinv, yerr = μ*sigma*pVinv*abs.(σ),
    linestyle="none", marker=".", elinewidth=0, capsize=3)
end

function plotmodelfield!(ax, iaxis, aem, ρ; color=nothing, alpha=1, model_lw=1, forward_lw=1)
    nfixed = aem.nfixed
    ax[iaxis].step(ρ, aem.z[nfixed+1:end], linewidth=model_lw, alpha=alpha)
    getfield!(ρ, aem)
    plotsoundingcurve(ax[iaxis+1], aem.F.dBzdt, aem.F.times; color, alpha, lw=forward_lw)
end    

function initmodelfield!(aem;  onesigma=true, figsize=(8,8))
    f, ax = plt.subplots(1, 2; figsize)
    if !isempty(aem.d)
        plotdata(ax[2], aem.d, aem.σ, aem.F.times; onesigma)
    end    
    ax
end    

function plotmodelfield!(aem, ρ; onesigma=true, color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,8), revax=true)
    plotmodelfield!(aem, [ρ]; onesigma, color, alpha, model_lw, forward_lw, figsize, revax) 
end  

function plotmodelfield!(aem, manyρ::Vector{T}; onesigma=true, 
        color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,8), revax=true) where T<:AbstractArray
    ax = initmodelfield!(aem; onesigma, figsize)
    for ρ in manyρ
        plotmodelfield!(ax, 1, aem, vec(ρ); alpha, model_lw, forward_lw, color)
    end
    ax[1].invert_yaxis()
    revax && ax[1].invert_xaxis()
    ax
end 

# noisy synthetic model making
function makenoisydata!(aem, ρ; 
        rseed=123, noisefrac=0.03, σ_halt=nothing, useML=false,
        onesigma=true, color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,8), revax=true)
    # σ_halt always assumed in Bfield units of pV
    getfield!(ρ, aem)
    f = aem.F.dBzdt
    σ_halt = isnothing(σ_halt) ? zeros(size(f)) : 1/pVinv*σ_halt/μ
    σ = sqrt.((noisefrac*abs.(f)).^2 + σ_halt.^2)
    Random.seed!(rseed)
    aem.d = f + σ.*randn(size(f))
    aem.σ = σ
    aem.useML = useML
    aem.ndata, aem.select = getndata(aem.d)
    aem.res, aem.J, aem.W = allocateJ(aem.F.dBzdt_J, aem.σ, aem.select, 
        aem.nfixed, length(aem.ρ), aem.F.calcjacobian)
    plotmodelfield!(aem, ρ; onesigma, color, alpha, model_lw, forward_lw, figsize, revax)
    nothing
end

function makeoperator(sounding::VTEMsoundingData; 
            useML         = false,
            zfixed        = [-1e5],
            ρfixed        = [1e12],
            extendfrac    = 1.06,
            dz            = 2.,
            ρbg           = 10, #linear ohm-m
            nlayers       = 40,
            ntimesperdecade = 10,
            nfreqsperdecade = 5,
            showgeomplot = false,
            calcjacobian = false,
            plotfield    = false
            )
    
    zall, znall, zboundaries = setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=showgeomplot)
    z, ρ, = makezρ(zboundaries; zfixed, ρfixed)
    aem = transD_GP.VTEM1DInversion.dBzdt(;d=sounding.data/μ, σ=sounding.noise/μ,
        times=sounding.times, ramp=sounding.ramp, ntimesperdecade, nfreqsperdecade,
        rTx=sounding.rTx, zTx=sounding.zTx, z, ρ, calcjacobian, useML, showgates=plotfield)
    
    plotfield && plotmodelfield!(aem, log10.(ρ[2:end]))
    aem, zall, znall
end

function plotwaveformgates(aem::dBzdt; figsize=(5,5))
    figure(;figsize)
    plot(aem.ramp[:,1]*1e6, aem.ramp[:,2], "-or")
    stem(aem.times*1e6, ones(length(aem.times)))
    ylabel("Amplitude")
    xlabel("time μs")
    title("Ramp and gates linear time")
end    

# for gradient based inversion
function loopacrosssoundings(soundings::Array{S, 1}, σstart, σ0; 
                            nsequentialiters   =-1,
                            zfixed             = [-1e5],
                            ρfixed             = [1e12],
                            zstart             = 0.0,
                            extendfrac         = 1.06,
                            dz                 = 2.,
                            ρbg                = 10,
                            nlayers            = 50,
                            ntimesperdecade    = 10,
                            nfreqsperdecade    = 5,
                            modelprimary       = false,
                            regtype            = :R1,
                            nstepsmax          = 10,
                            ntries             = 6,
                            target             = nothing,
                            lo                 = -3.,
                            hi                 = 1.,
                            λ²min              = 0,
                            λ²max              = 8,
                            λ²frac             = 4,
                            β²                 = 0.,
                            ntestdivsλ²        = 50,
                            αmin               = -4, 
                            αmax               = 0, 
                            αfrac              = 4, 
                            ntestdivsα         = 32,
                            regularizeupdate   = false,
                            knownvalue         = 0.7,
                            firstvalue         = :last,
                            κ                  = GP.Mat52(),
                            breakonknown       = true,
                            dobo               = false,
                            compresssoundings  = true,
                            zipsaveprefix      = "",        
                            ) where S<:Sounding

    @assert nsequentialiters  != -1
    nparallelsoundings = nworkers()
    nsoundings = length(soundings)
    zall, = setupz(zstart, extendfrac, dz=dz, n=nlayers) # needed for sounding compression
    for iter = 1:nsequentialiters
        if iter<nsequentialiters
            ss = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
        else
            ss = (iter-1)*nparallelsoundings+1:nsoundings
        end
        @info "soundings in loop $iter of $nsequentialiters", ss
        pids = workers()
        @sync for (i, s) in enumerate(ss)
            aem, = makeoperator(    soundings[s],
                                    zfixed = zfixed,
                                    ρfixed = ρfixed,
                                    zstart = zstart,
                                    extendfrac = extendfrac,
                                    dz = dz,
                                    ρbg = ρbg,
                                    nlayers = nlayers,
                                    modelprimary = modelprimary,
                                    ntimesperdecade = ntimesperdecade,
                                    nfreqsperdecade = nfreqsperdecade,
                                    calcjacobian = true)

            fname = soundings[s].sounding_string*"_gradientinv.dat"
            σstart_, σ0_ = map(x->x*ones(length(aem.ρ)-1), [σstart, σ0])
            @async remotecall_wait(gradientinv, pids[i], σstart_, σ0_, aem,
                                                regtype            = regtype         ,              
                                                nstepsmax          = nstepsmax       ,              
                                                ntries             = ntries          ,              
                                                target             = target          ,              
                                                lo                 = lo              ,              
                                                hi                 = hi              ,              
                                                λ²min              = λ²min           ,              
                                                λ²max              = λ²max           ,              
                                                λ²frac             = λ²frac          ,              
                                                ntestdivsλ²        = ntestdivsλ²     ,              
                                                αmin               = αmin            ,              
                                                αmax               = αmax            ,              
                                                αfrac              = αfrac           ,
                                                β²                 = β²              ,
                                                ntestdivsα         = ntestdivsα      ,              
                                                regularizeupdate   = regularizeupdate,              
                                                knownvalue         = knownvalue      ,              
                                                firstvalue         = firstvalue      ,              
                                                κ                  = κ               ,              
                                                breakonknown       = breakonknown    ,              
                                                dobo               = dobo            ,
                                                fname              = fname           ) 
                

        end # @sync
        isfirstparalleliteration = iter == 1 ? true : false
        compresssoundings && compress(soundings[ss[1]:ss[end]], zall, 
            isfirstparalleliteration = isfirstparalleliteration, prefix=zipsaveprefix)
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end

end