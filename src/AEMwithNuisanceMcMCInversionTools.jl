module AEMwithNuisanceMcMCInversionTools
using ..AbstractOperator, ..CommonToAll
import ..AbstractOperator.makeoperator
import ..AbstractOperator.loopacrossAEMsoundings
import ..AbstractOperator.plotmodelfield!
import ..AbstractOperator.makebounds
import ..AbstractOperator.getoptnfromexisting
import ..AbstractOperator.getnufromsounding

import ..Options, ..OptionsStat, ..OptionsNuisance
export makeAEMoperatorandnuisanceoptions, loopacrossAEMsoundings, summaryAEMnuisanceimages
import ..main # McMC function
using ..SoundingDistributor
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

# Driver code for McMC inversion with no nuisances
function loopacrossAEMsoundings(soundings::Array{S, 1}, aem_in::Operator1D, opt_in::Options, optn_in::OptionsNuisance;
                            Tmax               = -1,
                            nsamples           = -1,
                            nchainsatone       =  1,
                            nchainspersounding = -1,
                            ppn                = -1) where S<:Sounding

    @assert ppn != -1
    @assert nchainspersounding != -1
    @assert nsamples != - 1
    @assert Tmax != -1

    nsoundings = length(soundings)
    nsequentialiters, nparallelsoundings = splittasks(soundings; nchainspersounding, ppn)
    opt = deepcopy(opt_in)
    
    for iter = 1:nsequentialiters
        ss = getss(iter, nsequentialiters, nparallelsoundings, nsoundings)
        @info "soundings in loop $iter of $nsequentialiters", ss
        @sync for (i, s) in Iterators.reverse(enumerate(ss))
            pids = getpids(i, nchainspersounding)
            @info "pids in sounding $s:", pids
            
            aem = makeoperator(aem_in, soundings[s])
            opt.fdataname = soundings[s].sounding_string*"_"
            optn = getoptnfromexisting(optn_in, opt, soundings[s])

            @async remotecall_wait(main, pids[1], opt, optn, aem, collect(pids[2:end]);
                                    Tmax, nsamples, nchainsatone)

        end # @sync
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end

# plotting stuff

function summaryAEMwithnuisanceimages(soundings::Array{S, 1}, opt_in::Options, optn_in::OptionsNuisance;
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
                        idx = nothing,
                        omitconvergence = false,
                        preferEright = false,
                        preferNright = false,
                        numsize=1,
                        labelnu = [""],
                        dpi = 300,
                        saveplot = true,
                        yl = nothing,
                        showplot = true,
                        showmean = false,
                        ) where S<: Sounding

    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)                   
    for i in 1:nlines
        a = linestartidx[i]
        b = i != nlines ?  linestartidx[i+1]-1 : length(soundings)
        summaryimages(soundings[a:b], opt_in, optn_in; qp1, qp2, burninfrac, zall,dz, dr, 
            fontsize, vmin, vmax, cmap, figsize, topowidth, idx=idx, omitconvergence, useML, 
            preferEright, showplot, preferNright, saveplot, yl, dpi, showmean, numsize, labelnu)
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
                        figsize=(6,10),
                        topowidth=2,
                        idx = nothing,
                        omitconvergence = false,
                        preferEright = false,
                        preferNright = false,
                        numsize=1,
                        labelnu = [""],
                        dpi = 300,
                        saveplot = true,
                        yl = nothing,
                        showplot = true,
                        showmean = false,
                        ) where S<: Sounding
    @assert !(preferNright && preferEright) # can't prefer both labels to the right
    pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd, zall, 
    nulow, numid, nuhigh, nunominal = summarypost(soundings, opt_in, optn_in;
                                        qp1, qp2, burninfrac, zall, useML)

    phgrid, plgrid, pmgrid, σmeangrid, ∇zmeangrid,
    ∇zsdgrid, gridx, gridz, topofine, R, Z = makesummarygrid(soundings, pl, pm, ph, ρmean,
                                            vdmean, vddev, zall, dz, dr=dr)

    lname = "Line $(soundings[1].linenum)"
    plotsummarygrids1(soundings, σmeangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname; qp1, qp2,
                        figsize, fontsize, cmap, vmin, vmax, 
                        topowidth, idx, omitconvergence, useML,
                        preferEright, preferNright, saveplot, showplot, dpi,
                        yl, showmean)  

    plotsummarygrids3(soundings, nuhigh, nulow, numid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname, nunominal, numsize, labelnu=labelnu, qp1=qp1, qp2=qp2,
        figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax, 
        topowidth=topowidth, idx=idx, useML=useML, yl=yl,
        preferEright=preferEright, preferNright=preferNright, saveplot=saveplot)

end

function summarypost(soundings::Vector{S}, opt_in::Options, optn_in::OptionsNuisance;
        qp1=0.05,
        qp2=0.95,
        burninfrac=0.5,
        zall = [-1.],
        useML = false) where S<:Sounding

    @assert length(zall) != 1

    linename = "_line_$(soundings[1].linenum)_summary.txt"
    fnames = ["rho_low", "rho_mid", "rho_hi", "rho_avg",
        "ddz_mean", "ddz_sdev", "phid_mean", "phid_sdev",
        "nu_low", "nu_mid", "nu_high"].*linename

    opt = deepcopy(opt_in)
    idxnotzero = optn_in.idxnotzero
    nunominal = zeros(length(idxnotzero), length(soundings))
    if isfile(fnames[1])
        @warn fnames[1]*" exists, reading stored values"
        pl, pm, ph, ρmean,
        vdmean, vddev, χ²mean, χ²sd, 
        nulow, numid, nuhigh = map(x->readdlm(x), fnames)
        for idx = 1:length(soundings)
            nunominal[:,idx] = getnufromsounding(soundings[idx])[idxnotzero]
        end                                 
    else
        # this is a dummy operator for plotting
        aem, = makeoperator(soundings[1])
        pl, pm, ph, ρmean, vdmean, vddev = map(x->zeros(length(zall), length(soundings)), 1:6)
        χ²mean, χ²sd = zeros(length(soundings)), zeros(length(soundings))
        nulow, numid, nuhigh  = map(x->zeros(length(idxnotzero), length(soundings)), 1:3)
        for idx = 1:length(soundings)
            opt.fdataname = soundings[idx].sounding_string*"_"
            optn = getoptnfromexisting(optn_in, opt, soundings[idx])
            @info "$idx out of $(length(soundings))\n"
            opt.xall[:] .= zall
            nunominal[:,idx] = mean(optn.bounds[idxnotzero,:], dims=2)
            pl[:,idx], pm[:,idx], ph[:,idx], ρmean[:,idx],
            vdmean[:,idx], vddev[:,idx] = CommonToAll.plot_posterior(aem, opt, burninfrac=burninfrac,
                                                    qp1=qp1, qp2=qp2,
                                                    doplot=false)
            h, nuquants = plot_posterior(aem, optn, burninfrac=burninfrac, doplot=false)
            nulow[:, idx]  .= nuquants[:, 1]
            numid[:, idx]  .= nuquants[:, 2]
            nuhigh[:, idx] .= nuquants[:, 3]
            χ² = 2*CommonToAll.assembleTat1(opt, :U, temperaturenum=1, burninfrac=burninfrac)
            ndata = useML ? sum(.!isnan.(soundings[idx].Hx_data)) : sum(.!isnan.(soundings[idx].Hx_data)) +
                    sum(.!isnan.(soundings[idx].Hz_data))
            χ²mean[idx] = mean(χ²)/ndata
            χ²sd[idx]   = std(χ²)/ndata
            if useML 
                χ²mean[idx] = exp(χ²mean[idx]-log(ndata))
                χ²sd[idx]   = exp(χ²sd[idx]-log(ndata)) # I think, need to check
            end    
        end
        # write in grid format
        for (fname, vals) in Dict(zip(fnames, [pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd, nulow, numid, nuhigh]))
            writedlm(fname, vals)
        end
        # write in x, y, z, rho format
        for (i, d) in enumerate([pl, pm, ph, ρmean])
            xyzrho = makearray(soundings, d, zall)
            writedlm(fnames[i][1:end-4]*"_xyzrho.txt", xyzrho)
        end
    end
    pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd, zall, nulow, numid, nuhigh, nunominal
end

function plotsummarygrids3(soundings, nuhigh, nulow, numid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname, nunominal, numsize=2; qp1=0.05, qp2=0.95,
    figsize=(10,10), fontsize=12, cmap="viridis", vmin=-2, vmax=0.5, 
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
    s[icol].plot(R, ones(length(R)), "--k")
    s[icol].fill_between(R, vec(χ²mean-χ²sd), vec(χ²mean+χ²sd), alpha=0.5)
    s[icol].set_ylabel(L"ϕ_d")
    titlestring = useML ? "Max likelihood variance adjustment" : "Data misfit"
    s[icol].set_title(titlestring)
    icol += 1
    
    icol = plotnuquant(nulow, numid, nuhigh, nunominal, s, R, icol, nrows, numsize, labelnu)
    
    # pmgrid repeated as mean not shown in this plot
    summaryconductivity(s, icol, f, soundings, pmgrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, ; qp1, qp2, fontsize, 
        cmap, vmin, vmax, topowidth, idx, omitconvergence=false, preferEright, preferNright, yl, showmean=false)

    saveplot && savefig(lname*"_with_nu.png", dpi=dpi)
end

function plotnuquant(nqlow, nqmid, nqhigh, nunominal, s, gridx, icol, nrows, ms=2, labelnu=[""])
    nnu = min(size(nqlow, 1), size(nunominal, 1)) # in case we've inverted a zero bounds nuisance by mistech...
    labelnu[1] == "" || @assert length(labelnu) == nnu
    for inu = 1:nnu
        s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1])
        s[icol].fill_between(gridx, nqlow[inu,:], nqhigh[inu,:], alpha=0.5)
        s[icol].plot(gridx, nqmid[inu,:])
        s[icol].plot(gridx, nunominal[inu,:], "o", markersize=ms)
        labelnu[1] == "" || s[icol].set_title(labelnu[inu])
        labelnu[1] == "" || s[icol].set_ylabel(labelnu[inu])
        icol += 1
    end
    icol    
end   

function plotindividualsoundings(soundings::Vector{S}, 
    aem_in::Operator1D, opt_in::Options, optn_in::OptionsNuisance, 
    idxplot::Vector{Int};
    zall = [-1.],
    burninfrac=0.5,
    nbins = 50,
    figsize  = (6,6),
    omittemp = true,
    showslope = false,
    plotmean = false,
    pdfclim = nothing,
    qp1=0.05,
    qp2=0.95,
    rseed = 123,
    computeforwards = false,
    nforwards = 50) where S<:Sounding

    @assert length(zall) != 1
    
    opt = deepcopy(opt_in)
    opt.xall[:] = zall
    for idx = 1:length(soundings)
        if in(idx, idxplot)
            @info "Sounding number: $idx"
            aem = makeoperator(aem_in, soundings[idx])
            opt.fdataname = soundings[idx].sounding_string*"_"
            optn = getoptnfromexisting(optn_in, opt, soundings[idx]) 
            getchi2forall(opt, alpha=0.8; omittemp) # chi2 errors
            CommonToAll.getstats(opt) # ARs for GP model
            CommonToAll.getstats(optn) # ARs for nuisances
            plot_posterior(aem, opt; burninfrac, nbins, figsize,
                            showslope, pdfclim, plotmean) # GP models
            ax = gcf().axes
            ax[1].invert_xaxis()
            plot_posterior(aem, optn; burninfrac, nbins, figsize) # nuisances
            if computeforwards
                m = assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)
                mn = CommonToAll.assemblenuisancesatT(optn, temperaturenum=1, burninfrac=burninfrac)
                Random.seed!(rseed)
                plotmodelfield!(aem, m[randperm(length(m))[1:nforwards]],
                                    mn[randperm(length(m))[1:nforwards],:])
            end
        end
    end
end

end # module

