module AEMwithNuisanceMcMCInversionTools
using ..AbstractOperator, ..CommonToAll
import ..AbstractOperator.makeoperator
import ..AbstractOperator.loopacrossAEMsoundings
import ..AbstractOperator.plotmodelfield!
import ..AbstractOperator.getndata
import ..AbstractOperator.makebounds
import ..AbstractOperator.getoptnfromexisting
import ..AbstractOperator.getnufromsounding
import ..AbstractOperator.summaryAEMimages
import ..AbstractOperator.plotindividualAEMsoundings

import ..Options, ..OptionsStat, ..OptionsNuisance
export makeAEMoperatorandnuisanceoptions

import ..main # McMC function
using ..SoundingDistributor
import ..DEBUGLEVEL_TDGP
using Distributed, Dates, Statistics, DelimitedFiles, PyPlot, Random

function make_tdgp_opt(sounding::Sounding;
                    rseed = nothing,
                    znall = [1],
                    fileprefix = "sounding",
                    nmin = 2,
                    nmax = 40,
                    K = GP.Mat32(),
                    demean = false,
                    sampledc = true,
                    sddc = 0.01,
                    sdpos = 0.05,
                    sdprop = 0.05,
                    fbounds = [-0.5 2.5],
                    λ = [2],
                    δ = 0.1,
                    pnorm = 2,
                    save_freq = 50,
                    nuisance_sdev   = [0.],
                    nuisance_bounds = [0. 0.],
                    updatenuisances = true,
                    C = nothing,
                    dispstatstoscreen = false,
                    restart = false
                    )
    sdev_dc = sddc*diff(fbounds, dims=2)[:]
    sdev_pos = [sdpos*abs(diff([extrema(znall)...])[1])]
    sdev_prop = sdprop*diff(fbounds, dims=2)[:]
    xall = permutedims(collect(znall))
    xbounds = permutedims([extrema(znall)...])

    history_mode = "w"
    restart && (history_mode = "a")

    updatenonstat = false
    needλ²fromlog = false
    if rseed != nothing
        Random.seed!(rseed)
    end

    opt = OptionsStat(fdataname = fileprefix*"_",
                            nmin = nmin,
                            nmax = nmax,
                            xbounds = xbounds,
                            fbounds = fbounds,
                            xall = xall,
                            λ = λ,
                            δ = δ,
                            demean = demean,
                            sampledc = sampledc,
                            sdev_dc = sdev_dc,
                            sdev_prop = sdev_prop,
                            sdev_pos = sdev_pos,
                            pnorm = pnorm,
                            quasimultid = false,
                            K = K,
                            save_freq = save_freq,
                            needλ²fromlog = needλ²fromlog,
                            updatenonstat = updatenonstat,
                            updatenuisances = updatenuisances,
                            dispstatstoscreen = dispstatstoscreen,
                            history_mode = history_mode
                            )

    bounds = makebounds(nuisance_bounds, sounding::Sounding)                   

    optn = OptionsNuisance(opt;
                    sdev = copy(nuisance_sdev),
                    bounds = bounds, C = C,
                    updatenuisances = updatenuisances)

                    opt, optn
                end
                
function makeAEMoperatorandnuisanceoptions(sounding::Sounding;
                                                    zfixed   = [-1e5],
                                                    ρfixed   = [1e12],
                                                    zstart = 0.0,
                                                    extendfrac = 1.06,
                                                    useML = false,
                                                    dz = 2.,
                                                    ρbg = 10,
                                                    nlayers = 40,
                                                    ntimesperdecade = 10,
                                                    nfreqsperdecade = 5,
                                                    showgeomplot = false,
                                                    plotfield = true,
                                                    nmin = 2,
                                                    nmax = 40,
                                                    K = GP.Mat32(),
                                                    demean = false,
                                                    sampledc = true,
                                                    sddc = 0.01,
                                                    sdpos = 0.05,
                                                    sdprop = 0.05,
                                                    fbounds = [-0.5 2.5],
                                                    λ = [2],
                                                    δ = 0.1,
                                                    save_freq = 50,
                                                    nuisance_sdev   = [0.],
                                                    nuisance_bounds = [0. 0.],
                                                    C = nothing,
                                                    updatenuisances = true,
                                                    dispstatstoscreen = false,
                                                    vectorsum = false,
                                                    restart = false
                                                    )
    
    aem, zall, znall, = makeoperator(sounding;
                                    zfixed, ρfixed, zstart, extendfrac,
                                    dz, ρbg, nlayers, ntimesperdecade,
                                    nfreqsperdecade, showgeomplot,
                                    plotfield, useML, vectorsum)
    
    opt, optn = make_tdgp_opt(sounding;
                                    znall = znall,
                                    fileprefix = sounding.sounding_string,
                                    nmin = nmin,
                                    nmax = nmax,
                                    K = K,
                                    demean = demean,
                                    sampledc = sampledc,
                                    sddc = sddc,
                                    sdpos = sdpos,
                                    sdprop = sdprop,
                                    fbounds = fbounds,
                                    save_freq = save_freq,
                                    λ = λ,
                                    δ = δ,
                                    nuisance_bounds = nuisance_bounds,
                                    nuisance_sdev = nuisance_sdev,
                                    updatenuisances = updatenuisances,
                                    C = C,
                                    dispstatstoscreen = dispstatstoscreen,
                                    restart
                                    )
    aem, opt, optn, zall
end

# Driver code for McMC inversion with nuisances
function loopacrossAEMsoundings(soundings::Array{S, 1}, aem_in::Operator1D, opt_in::Options, optn_in::OptionsNuisance;
                            Tmax               = -1,
                            nsamples           = -1,
                            nchainsatone       =  1,
                            nchainspersounding = -1,
                            ppn                = -1,
                            nominaltime        = nothing) where S<:Sounding

    @assert ppn != -1
    @assert nchainspersounding != -1
    @assert nsamples != - 1
    @assert Tmax != -1

    nsoundings = length(soundings)
    nsequentialiters, nparallelsoundings = splittasks(soundings; nchainspersounding, ppn)
    
    for iter = 1:nsequentialiters
        ss = getss(iter, nsequentialiters, nparallelsoundings, nsoundings)
        @info "soundings in loop $iter of $nsequentialiters", ss
        t2 = time()
        @sync for (i, s) in Iterators.reverse(enumerate(ss))
            pids = getpids(i, nchainspersounding)
            (DEBUGLEVEL_TDGP > 0) && @info("pids in sounding $s are $pids")
            aem = makeoperator(aem_in, soundings[s])
            opt = deepcopy(opt_in) # this is a big bugfix ... https://discourse.julialang.org/t/unexpected-behavior-of-async/91921
            opt.fdataname = soundings[s].sounding_string*"_"
            optn = getoptnfromexisting(optn_in, opt, soundings[s])

            @async remotecall_wait(main, pids[1], opt, optn, aem, collect(pids[2:end]);
                                    Tmax, nsamples, nchainsatone, nominaltime)

        end # @sync
        dt = time() - t2 #seconds
        t2 = time()
        @info "done $iter out of $nsequentialiters at $(Dates.now()) in $dt sec"
    end
end

# plotting stuff

function summaryAEMimages(soundings::Array{S, 1}, opt_in::Options, optn_in::OptionsNuisance;
                        zall=[-1.],
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        useML = false,
                        dz = 2.,
                        dr = 10,
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap="turbo",
                        figsize=(6,10),
                        topowidth=2,
                        lnames = [], # array of lines
                        idx = [], # array of arrrays per line
                        omitconvergence = false,
                        preferEright = false,
                        preferNright = false,
                        numsize=1,
                        labelnu = [""],
                        dpi = 300,
                        saveplot = true,
                        yl = nothing,
                        showplot = true,
                        logscale = true,
                        showmean = false,
                        vectorsum = false,
                        Rmax = nothing,
                        ) where S<: Sounding
    compatidxwarn(idx, lnames)
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)                   
    for i in 1:nlines
        a, b = linestartend(linestartidx, i, nlines, soundings)
        continueflag, idspec = docontinue(lnames, idx, soundings, a, b)
        continueflag && continue
        summaryimages(soundings[a:b], opt_in, optn_in; qp1, qp2, burninfrac, zall,dz, dr, vectorsum,
            fontsize, vmin, vmax, cmap, figsize, topowidth, idx=idspec, omitconvergence, useML, logscale,
            preferEright, showplot, preferNright, saveplot, yl, dpi, showmean, numsize, labelnu, Rmax)
    end
    nothing  
end                        

function summaryimages(soundings::Array{S, 1}, opt_in::Options, optn_in::OptionsNuisance;
                        zall=[-1.],
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        useML = false,
                        dz = 2.,
                        dr = 10,
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap="turbo",
                        figsize = (6,10),
                        bigfigsize = figsize,
                        topowidth = 2,
                        idx = [],
                        omitconvergence = false,
                        preferEright = false,
                        preferNright = false,
                        numsize=1,
                        labelnu = [""],
                        logscale = true,
                        dpi = 300,
                        saveplot = true,
                        yl = nothing,
                        showplot = true,
                        showmean = false,
                        vectorsum = false,
                        Rmax = nothing,
                        ) where S<: Sounding
    @assert !(preferNright && preferEright) # can't prefer both labels to the right
    pl, pm, ph, ρmean,  χ²mean, χ²sd,  
    nulow, numid, nuhigh, nunominal = summarypost(soundings, opt_in, optn_in; vectorsum,
                                        qp1, qp2, burninfrac, zall, useML)

    phgrid, plgrid, pmgrid, σmeangrid, 
    gridx, gridz, topofine, R, Z = makesummarygrid(soundings, pl, pm, ph, ρmean,
                                            zall, dz; dr)

    lname = "Line_$(soundings[1].linenum)"
    plotsummarygrids1(soundings, σmeangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname; qp1, qp2,
                        figsize, fontsize, cmap, vmin, vmax, 
                        topowidth, idx, omitconvergence, useML, Rmax, logscale,
                        preferEright, preferNright, saveplot, showplot, dpi,
                        yl, showmean)  

    plotsummarygrids3(soundings, nuhigh, nulow, numid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname, nunominal, numsize; labelnu=labelnu, qp1=qp1, qp2=qp2,
        figsize=bigfigsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax, 
        topowidth=topowidth, idx=idx, useML=useML, yl=yl, logscale=logscale, showplot,
        preferEright=preferEright, preferNright=preferNright, saveplot=saveplot)

end

function summarypost(soundings::Vector{S}, opt_in::Options, optn_in::OptionsNuisance;
        qp1=0.05,
        qp2=0.95,
        burninfrac=0.5,
        zall = [-1.],
        vectorsum = false,
        useML = false) where S<:Sounding

    @assert length(zall) != 1

    linename = "_line_$(soundings[1].linenum)_summary.txt"
    fnames = ["rho_low", "rho_mid", "rho_hi", "rho_avg",
        "phid_mean", "phid_sdev",
        "nu_low", "nu_mid", "nu_high"].*linename
    idxnotzero = optn_in.idxnotzero
    # this is a debug for unfinished soundings
    a = Vector{Any}(undef, 10)
    for i in 1:4
        a[i] = -999*ones(length(zall))
    end
    a[5] = -999.
    a[6] = -999.
    for i in 7:10
        a[i] = -999*ones(length(idxnotzero))
    end
    if isfile(fnames[1])
        nunominal = zeros(length(idxnotzero), length(soundings))
        @warn fnames[1]*" exists, reading stored values"
        pl, pm, ph, ρmean,
        χ²mean, χ²sd, 
        nulow, numid, nuhigh = map(x->readdlm(x), fnames)
        for idx = 1:length(soundings)
            nunominal[:,idx] = getnufromsounding(soundings[idx])[idxnotzero]
        end                                 
    else
        outputs = reduce(hcat, pmap(sounding->processonesounding(opt_in, 
                                optn_in, sounding, zall, burninfrac, qp1, qp2, idxnotzero, useML, vectorsum), soundings;
                                on_error=ex->a))
        pl, pm, ph, ρmean,  
        χ²mean, χ²sd, nulow, 
        numid, nuhigh, nunominal = map(x->reduce(hcat, x), eachrow(outputs))
        χ²mean, χ²sd = map(x->vec(x), (χ²mean, χ²sd))
        # write in grid format
        for (fname, vals) in Dict(zip(fnames, [pl, pm, ph, ρmean, χ²mean, χ²sd, nulow, numid, nuhigh]))
            writedlm(fname, vals)
        end
        # write in x, y, z, rho format
        for (i, d) in enumerate([pl, pm, ph, ρmean])
            xyzrho = makearray(soundings, d, zall)
            writedlm(fnames[i][1:end-4]*"_xyzrho.txt", xyzrho)
        end
    end
    pl, pm, ph, ρmean, χ²mean, χ²sd, nulow, numid, nuhigh, nunominal
end

function processonesounding(opt_in::Options, optn_in::OptionsNuisance, sounding::Sounding, zall, burninfrac, qp1, qp2, idxnotzero, useML, vectorsum)
    opt = deepcopy(opt_in)
    opt.fdataname = sounding.sounding_string*"_"
    optn = getoptnfromexisting(optn_in, opt, sounding)
    @info ("processing FID $(sounding.fid)")
    opt.xall[:] .= zall
    nunominal = mean(optn.bounds[idxnotzero,:], dims=2)
    aem, = makeoperator(sounding) # this is a dummy operator for plotting
    pl, pm, ph, ρmean, = CommonToAll.plot_posterior(aem, opt; burninfrac, qp1, qp2, doplot=false)
    _, nuquants = plot_posterior(aem, optn; burninfrac, doplot=false)
    nulow, numid, nuhigh = nuquants[:,1], nuquants[:,2], nuquants[:,3]
    χ² = 2*CommonToAll.assembleTat1(opt, :U; temperaturenum=1, burninfrac)
    ndata = getndata(sounding, vectorsum)
    χ²mean = mean(χ²)/ndata
    χ²sd   = std(χ²)/ndata
    if sounding.forceML
        χ²mean = 1.
        χ²sd   = 0.
    elseif useML
        # same ML factor for Hx and Hz beware
        χ²mean = exp(χ²mean-log(ndata))
        χ²sd   = exp(χ²sd-log(ndata)) # I think, need to check
    end
    [pl, pm, ph, ρmean, χ²mean, χ²sd, nulow, numid, nuhigh, nunominal]
end

function plotsummarygrids3(soundings, nuhigh, nulow, numid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname, nunominal, numsize=2; qp1=0.05, qp2=0.95,
    figsize=(10,10), fontsize=12, cmap="viridis", vmin=-2, vmax=0.5, logscale=true, showplot = true,
    topowidth=2, idx=nothing, useML=false, preferEright=false, preferNright=false, labelnu=[""], yl = nothing, dpi=300,
    saveplot=true)
    
    # get the finely interpolated nuisances
    # nuhighfine, nulowfine, numidfine = map(X->gridpoints(R, gridx, X), [nuhigh, nulow, numid])
    
    dr = diff(gridx)[1]
    nnu = min(size(nulow, 1), size(nunominal, 1)) # in case we've inverted a zero bounds nuisance by mistech...
    nrows = 1 + nnu + 3 + 1 # add the number of nuisances == no. of rows in nuhigh and 1 for chi2 and 1 for colorbar, not showing mean
    height_ratios = [0.4ones(1+nnu)...,1,1,1,0.1]
    f, s = plt.subplots(nrows, 1, gridspec_kw=Dict("height_ratios" => height_ratios),
    figsize=figsize)
    f.suptitle(lname*" Δx=$dr m, Fids: $(length(R))", fontsize=fontsize)
    icol = 1
    s[icol].plot(R, χ²mean)
    linecolor = [s[icol].get_facecolor()...]
    linecolor[1:3] = abs.(linecolor[1:3] .- 1)
    s[icol].plot(R, ones(length(R)), "--", color=linecolor)
    s[icol].fill_between(R, vec(χ²mean-χ²sd), vec(χ²mean+χ²sd), alpha=0.5)
    s[icol].set_ylabel(L"ϕ_d")
    titlestring = useML ? "Max likelihood variance adjustment" : "Data misfit"
    s[icol].set_title(titlestring)
    logscale && s[icol].set_yscale("log")
    icol += 1
    
    icol = plotnuquant(nulow, numid, nuhigh, nunominal, s, R, icol, nrows, numsize, labelnu)
    
    # pmgrid repeated as mean not shown in this plot
    summaryconductivity(s, icol, f, soundings, pmgrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, ; qp1, qp2, fontsize, 
        cmap, vmin, vmax, topowidth, idx, omitconvergence=false, preferEright, preferNright, yl, showmean=false)

    saveplot && savefig(lname*"_with_nu.png", dpi=dpi)
    showplot || close(f)
end


function plotnuquant(nqlow, nqmid, nqhigh, nunominal, s, gridx, icol, nrows, ms=2, labelnu=[""])
    nnu = min(size(nqlow, 1), size(nunominal, 1)) # in case we've inverted a zero bounds nuisance by mistech...
    labelnu[1] == "" || @assert length(labelnu) == nnu
    for inu = 1:nnu
        s[icol].sharex(s[icol-1])
        s[icol].fill_between(gridx, nqlow[inu,:], nqhigh[inu,:], alpha=0.5)
        s[icol].plot(gridx, nqmid[inu,:])
        s[icol].plot(gridx, nunominal[inu,:], "o", markersize=ms)
        labelnu[1] == "" || s[icol].set_title(labelnu[inu])
        labelnu[1] == "" || s[icol].set_ylabel(labelnu[inu])
        icol += 1
    end
    icol    
end   

function plotindividualAEMsoundings(soundings::Vector{S}, 
    aem_in::Operator1D, opt_in::Options, optn_in::OptionsNuisance, 
    idxplot;
    zall = [-1.],
    burninfrac=0.5,
    nbins = 50,
    figsize  = (6,6),
    omittemp = true,
    showslope = false,
    plotmean = false,
    pdfclim = nothing,
    model_lw = 1, 
    forward_lw = 1,
    qp1=0.05,
    qp2=0.95,
    linecolor = nothing,
    alpha = 1.,
    rseed = 123,
    lnames = [],
    usekde = false,
    computeforwards = false,
    nforwards = 20) where S<:Sounding

    compatidxwarn(idxplot, lnames)
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)        
    @assert length(zall) != 1
    opt = deepcopy(opt_in)
    opt.xall[:] = zall
    for i in 1:nlines
        a, b = linestartend(linestartidx, i, nlines, soundings)
        continueflag, idspec = docontinue(lnames, idxplot, soundings, a, b)
        continueflag && continue
        for idx in idspec
            @info "Sounding number: $idx"
            aem = makeoperator(aem_in, soundings[a:b][idx])
            opt.fdataname = soundings[a:b][idx].sounding_string*"_"
            optn = getoptnfromexisting(optn_in, opt, soundings[a:b][idx]) 
            getchi2forall(opt, alpha=0.8; omittemp) # chi2 errors
            CommonToAll.getstats(opt) # ARs for GP model
            gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
            nicenup(gcf())
            CommonToAll.getstats(optn) # ARs for nuisances
            gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
            nicenup(gcf())
            plot_posterior(aem, opt; burninfrac, nbins, figsize, qp1, qp2,
                            showslope, pdfclim, plotmean, usekde) # GP models
            ax = gcf().axes
            gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
            ax[1].invert_xaxis()
            nicenup(gcf())
            plot_posterior(aem, optn; burninfrac, nbins, figsize, qp1, qp2) # nuisances
            gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
            nicenup(gcf())
            if computeforwards
                m = assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)
                mn = CommonToAll.assemblenuisancesatT(optn, temperaturenum=1, burninfrac=burninfrac)
                Random.seed!(rseed)
                randidx = randperm(length(m))
                plotmodelfield!(aem, m[randidx[1:nforwards]], mn[randidx[1:nforwards],:]; model_lw, forward_lw, color=linecolor, alpha)
                gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
                nicenup(gcf())
            end
        end
    end    
end

# legacy backwards compat
plotindividualAEMsoundingswithnuisance = plotindividualAEMsoundings
plotindividualsoundings = plotindividualAEMsoundings
summaryAEMwithnuisanceimages = summaryAEMimages
end # module

