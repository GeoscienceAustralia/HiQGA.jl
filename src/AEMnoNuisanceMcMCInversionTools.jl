module AEMnoNuisanceMcMCInversionTools
using ..AbstractOperator, ..CommonToAll
import ..AbstractOperator.makeoperator
import ..AbstractOperator.loopacrossAEMsoundings
import ..AbstractOperator.plotmodelfield!
import ..AbstractOperator.getndata
import ..Options, ..OptionsStat
import ..AbstractOperator.summaryAEMimages
import ..AbstractOperator.plotindividualAEMsoundings

export makeAEMoperatorandoptions, loopacrossAEMsoundings, summaryAEMimages, plotindividualAEMsoundings
import ..main # McMC function
import ..DEBUGLEVEL_TDGP
using ..SoundingDistributor
using Distributed, Dates, Statistics, DelimitedFiles, PyPlot, Random
# plotting stuff
function summaryAEMimages(soundings::Array{S, 1}, opt::Options;
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        dz = -1,
                        dr = 10,
                        zall=[-1.],
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap="turbo",
                        figsize=(6,10),
                        topowidth=2,
                        lnames = [], # array of lines
                        idx = [], # array of arrrays per line
                        omitconvergence = false,
                        useML = false,
                        preferEright = false,
                        preferNright = false,
                        saveplot = true,
                        yl = nothing,
                        showplot = true,
                        showmean = false,
                        logscale = true,
                        dpi = 300) where S<:Sounding
    compatidxwarn(idx, lnames)
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)                   
    for i in 1:nlines
        a, b = linestartend(linestartidx, i, nlines, soundings)
        continueflag, idspec = docontinue(lnames, idx, soundings, a, b)
        continueflag && continue
        summaryimages(soundings[a:b], opt; qp1, qp2, burninfrac, zall,dz, dr, 
            fontsize, vmin, vmax, cmap, figsize, topowidth, idx=idspec, omitconvergence, useML, 
            preferEright, showplot, preferNright, saveplot, yl, dpi, showmean, logscale)
    end
    nothing    
end

function summaryimages(soundings::Array{S, 1}, opt::Options;
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        dz = -1,
                        dr = 10,
                        zall = [-1.],
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap ="turbo",
                        figsize = (6,10),
                        topowidth=2,
                        idx = [],
                        omitconvergence = false,
                        useML = false,
                        preferEright = false,
                        preferNright = false,
                        saveplot = false,
                        yl = nothing,
                        logscale = true,
                        showplot = true,
                        showmean = false,
                        dpi = 300) where S<:Sounding
    @assert !(preferNright && preferEright) # can't prefer both labels to the right
    pl, pm, ph, ρmean, χ²mean, χ²sd  = summarypost(soundings, opt; zall, qp1, qp2, burninfrac, useML)

    phgrid, plgrid, pmgrid, σmeangrid,
    gridx, gridz, topofine, R, Z = makesummarygrid(soundings, pl, pm, ph, ρmean,
                                                        zall, dz; dr)

    lname = "Line_$(soundings[1].linenum)"
    plotsummarygrids1(soundings, σmeangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname; qp1, qp2,
                        figsize, fontsize, cmap, vmin, vmax, 
                        topowidth, idx, omitconvergence, useML,
                        preferEright, preferNright, saveplot, showplot, dpi,
                        yl, showmean, logscale)                  
end

function summarypost(soundings::Vector{S}, opt::Options;
            qp1=0.05,
            qp2=0.95,
            burninfrac=0.5,
            zall = [-1.],
            useML=false) where S<:Sounding

    @assert length(zall) != 1

    linename = "_line_$(soundings[1].linenum)_summary.txt"
    fnames = ["rho_low", "rho_mid", "rho_hi", "rho_avg",
              "phid_mean", "phid_sdev"].*linename
    # this is a debug for unfinished soundings
    a = Vector{Any}(undef, 6)
    for i in 1:4
        a[i] = -999*ones(length(zall))
    end
    a[5] = -999.
    a[6] = -999.
    if isfile(fnames[1])
        @warn fnames[1]*" exists, reading stored values"
        pl, pm, ph, ρmean,
        χ²mean, χ²sd = map(x->readdlm(x), fnames)
    else
        outputs = reduce(hcat, pmap(sounding->processonesounding(opt, 
                                sounding, zall, burninfrac, qp1, qp2, useML), soundings))
        pl, pm, ph, ρmean,  
        χ²mean, χ²sd = map(x->reduce(hcat, x), eachrow(outputs))
        χ²mean, χ²sd = map(x->vec(x), (χ²mean, χ²sd))
        # write in grid format
        for (fname, vals) in Dict(zip(fnames, [pl, pm, ph, ρmean, χ²mean, χ²sd]))
            writedlm(fname, vals)
        end
        # write in x, y, z, rho format
        for (i, d) in enumerate([pl, pm, ph, ρmean])
            xyzrho = makearray(soundings, d, zall)
            writedlm(fnames[i][1:end-4]*"_xyzrho.txt", xyzrho)
        end
    end
    pl, pm, ph, ρmean, χ²mean, χ²sd
end

function processonesounding(opt_in::Options, sounding::Sounding, zall, burninfrac, qp1, qp2, useML)
    opt = deepcopy(opt_in)
    opt.fdataname = sounding.sounding_string*"_"
    @info ("processing FID $(sounding.fid)")
    opt.xall[:] .= zall
    aem, = makeoperator(sounding) # this is a dummy operator for plotting
    pl, pm, ph, ρmean, = CommonToAll.plot_posterior(aem, opt; burninfrac, qp1, qp2, doplot=false)
    χ² = 2*CommonToAll.assembleTat1(opt, :U, temperaturenum=1, burninfrac=burninfrac)
    ndata = getndata(sounding)
    χ²mean = mean(χ²)/ndata
    χ²sd   = std(χ²)/ndata
    if hasproperty(sounding, :forceML) 
        if sounding.forceML
            χ²mean = 1.
            χ²sd   = 0.
        end     
    elseif useML
        # this is approximate as HM and LM have different ML factors sampled 
        χ²mean = exp(χ²mean-log(ndata))
        χ²sd   = exp(χ²sd-log(ndata)) # I think, need to check
    end
    [pl, pm, ph, ρmean, χ²mean, χ²sd]
end

function getndata(d)
    select = .!isnan.(d)
    ndata  = sum(select)
    ndata, select
end

function plotindividualAEMsoundings(soundings::Vector{S}, aem_in::Operator1D, opt_in::Options, idxplot; # idxplot is an array of arrays
    lnames = [], # array of lines
    zall=[-1.],
    burninfrac=0.5,
    nbins = 50,
    figsize  = (6,6),
    omittemp = true,
    showslope = false,
    plotmean = false,
    pdfclim = nothing,
    qp1 = 0.05,
    qp2 = 0.95,
    model_lw = 1, 
    forward_lw = 1,
    linecolor = nothing,
    alpha = 1.,
    rseed = 123,
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
            getchi2forall(opt, alpha=0.8, omittemp=omittemp)
            gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
            nicenup(gcf())
            CommonToAll.getstats(opt)
            gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
            nicenup(gcf())
            plot_posterior(aem, opt; burninfrac, nbins, figsize, 
                    showslope, pdfclim, plotmean, qp1, qp2, usekde)
            ax = gcf().axes
            gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
            ax[1].invert_xaxis()
            nicenup(gcf())
            if computeforwards
                M = assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)
                Random.seed!(rseed)
                plotmodelfield!(aem, M[randperm(length(M))[1:nforwards]]; model_lw, forward_lw, color=linecolor, alpha)
                gcf().suptitle("Line $(soundings[a].linenum) index:$idx")
                nicenup(gcf())
            end            
        end
    end    
end

# stuff needed for McMC driver code
function make_tdgp_opt(;
                    rseed = nothing,
                    znall = znall,
                    fileprefix = "sounding",
                    nmin = 2,
                    nmax = 40,
                    K = GP.OrstUhn(),
                    demean = true,
                    sdpos = 0.05,
                    sdprop = 0.05,
                    sddc = 0.008,
                    sampledc = false,
                    fbounds = [-0.5 2.5],
                    λ = [2],
                    δ = 0.1,
                    pnorm = 2,
                    save_freq = 25,
                    restart = false
                    )
    sdev_pos = [sdpos*abs(diff([extrema(znall)...])[1])]
    sdev_prop = sdprop*diff(fbounds, dims=2)[:]
    sdev_dc = sddc*diff(fbounds, dims=2)[:]
    xall = permutedims(collect(znall))
    xbounds = permutedims([extrema(znall)...])

    history_mode = "w"
	restart && (history_mode = "a")

    if rseed != nothing
        Random.seed!(rseed)
    end

    opt = OptionsStat(;fdataname = fileprefix*"_",
        nmin, nmax, xbounds, fbounds, xall, λ, δ,
        demean, sdev_prop, sdev_pos, sdev_dc,
        sampledc, pnorm, quasimultid = false, K,
        save_freq, needλ²fromlog=false, updatenonstat=false,
        dispstatstoscreen = false, history_mode)
    opt
end

function makeAEMoperatorandoptions(sounding::Sounding;
                        nmin = 2,
                        nmax = 40,
                        K = GP.OrstUhn(),
                        demean = false,
                        sdpos = 0.05,
                        sdprop = 0.05,
                        sddc = 0.008,
                        sampledc = true,
                        fbounds = [-0.5 2.5],
                        λ = [2],
                        δ = 0.1,
                        save_freq = 25,
                        restart = false,
                        zfixed   = [-1e5],
                        ρfixed   = [1e12],
                        zstart = 0.0,
                        extendfrac = 1.06,
                        dz = 2.,
                        ρbg = 10,
                        nlayers = 40,
                        ntimesperdecade = 10,
                        nfreqsperdecade = 5,
                        showgeomplot = false,
                        plotfield = false,
                        useML = false,
                        modelprimary = false
                        )
    aem, zall, znall,  = makeoperator(sounding;
                        zfixed, ρfixed, zstart, extendfrac,
                        dz = dz, ρbg, nlayers, ntimesperdecade,
                        nfreqsperdecade, useML, showgeomplot,
                        modelprimary, plotfield)

    opt = make_tdgp_opt(znall = znall,
                        fileprefix = sounding.sounding_string,
                        nmin = nmin,
                        nmax = nmax,
                        K = K,
                        demean = demean,
                        sdpos = sdpos,
                        sdprop = sdprop,
                        sddc = sddc,
                        sampledc = sampledc,
                        fbounds = fbounds,
                        save_freq = save_freq,
                        λ = λ,
                        δ = δ,
                        restart = restart
                        )
    aem, opt, zall
end

# Driver code for McMC inversion with no nuisances
function loopacrossAEMsoundings(soundings::Array{S, 1}, aem_in::Operator1D, opt_in::Options;
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
            opt = deepcopy(opt_in)
            opt.fdataname = soundings[s].sounding_string*"_"

            @async remotecall_wait(main, pids[1], opt, aem, collect(pids[2:end]);
                                    Tmax, nsamples, nchainsatone, nominaltime)

        end # @sync
        dt = time() - t2 #seconds
        t2 = time()
        @info "done $iter out of $nsequentialiters at $(Dates.now()) in $dt sec"
    end
end

end