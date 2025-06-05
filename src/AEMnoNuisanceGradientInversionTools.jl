module AEMnoNuisanceGradientInversionTools
using Distributed, Dates, Printf, PyPlot, DelimitedFiles, StatsBase, Random
using ..AbstractOperator, ..CommonToAll, ..GP, ..SoundingDistributor
import ..AbstractOperator.makeoperator
import ..AbstractOperator.Sounding
import ..AbstractOperator.returnforwrite
import ..AbstractOperator.loopacrossAEMsoundings
import ..AbstractOperator.plotconvandlast
import ..AbstractOperator.plotmodelfield!
import ..AbstractOperator.getresidual
import ..gradientinv
export plotconvandlast, loopacrossAEMsoundings
# for deterministic inversions, read in
function compress(soundings, zall; prefix="", rmfile=true, isfirstparalleliteration=false)
    fname = soundings[1].sounding_string*"_gradientinv.dat"    
    !isfile(fname) && throw(AssertionError("file $fname does not exist perhaps soundings already zipped?"))
    fout = prefix == "" ? "zipped.dat" : prefix*"_zipped.dat"
    isfirstparalleliteration && isfile(fout) && throw(AssertionError("Zipped file "*fout*" exists, will not overwrite!"))
    for (i, s) in enumerate(soundings)
        fname = s.sounding_string*"_gradientinv.dat"
        A = readdlm(fname)
        ϕd = A[end,2]
        σgrid = vec(A[end,3:end])
        elinonerow = [returnforwrite(s)..., vec(zall), σgrid, ϕd]
        nelinonerow = length(elinonerow)
        writenames = [string.(s.writefields)..., "zcenter", "log10_cond", "phid_err"]
        sfmt = fill("%15.3f", nelinonerow)
        ϕd > 1e4 && (ϕd = 1e4 )
        iomode = (isfirstparalleliteration && i==1) ? "w" : "a"
        writeasegdat(elinonerow, sfmt, fout[1:end-4], iomode)
        rmfile && rm(fname)
        if isfirstparalleliteration && i == 1
            channel_names = [writenames, fill("", nelinonerow), writenames]
            writeasegdfnfromonerow(elinonerow, channel_names, sfmt, fout[1:end-4])
            dfn2hdr(fout[1:end-4]*".dfn")
        end
    end
end

# plot the convergence and the result
function plotconvandlast(soundings, delr, delz; 
        zall=[-1.], cmapσ="turbo", vmin=-2.5, vmax=0.5, fontsize=12,
        figsize=(20,5),
        topowidth=1,
        preferEright = false,
        preferNright = false,
        saveplot = true,
        showplot = true,
        postfix = "",
        prefix = "",
        markersize = 2,
        logscale = false,
        lnames = [],
        idx = [], # array of arrrays per line
        yl = nothing,
        plotforward = false,
        aem_in = nothing,
        nnu = 0,
        dophiplot = false,
        doreshist = false,
        calcresiduals = false,
        dpi=400)
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)
    nlayers = length(zall)
    fnamecheck = soundings[1].sounding_string*"_gradientinv.dat"
    isfile(fnamecheck) && compress(soundings, zall) # write everything in one file if not done yet
    fzipped = prefix == "" ? "zipped.dat" : prefix*"_zipped.dat"
    A = readdlm(fzipped)
    σ = A[:,end-nlayers-nnu:end-nnu-1] # so we can plot TEMPEST and SPECTREM similar to heli
    ϕd = A[:,end]
    res = calcresiduals ? Vector{Array{Float64, 2}}(undef, nlines) : fill(zeros(0), nlines)
    for i in 1:nlines
        a, b = linestartend(linestartidx, i, nlines, soundings)
        continueflag, idspec = docontinue(lnames, idx, soundings, a, b)
        continueflag && continue
        if calcresiduals
            nsoundings = b-a+1
            res[i] = zeros(length(aem_in.res), nsoundings) # won't work if different number of data per sounding
            for (id, s) in enumerate(soundings[a:b])
                aem = makeoperator(aem_in, s)
                m = vec(σ[a:b,:][id,:]) # log 10 σ
                getresidual(aem, m) # used for grad inversions so log10 cond
                res[i][:,id] .= aem.W*aem.res
            end
        end
        if plotforward && !isempty(idspec) && !isnothing(aem_in)
            for id in idspec
                @info "ϕd is $(ϕd[a:b][id])"
                aem = makeoperator(aem_in, soundings[a:b][id])
                m = -vec(σ[a:b,:][id,:]) #log 10 ρ
                plotmodelfield!(aem, m)
                gcf().suptitle("Line $(soundings[a].linenum) index:$id")
                nicenup(gcf())
            end    
        end   
        # skip for debugging and highly decimated surveys
        a == b && continue
        plotconvandlasteachline(soundings[a:b], view(σ, a:b, :)', view(ϕd, a:b), delr, delz, res[i]; 
            zall = zall, idx=idspec, yl=yl,
            cmapσ=cmapσ, vmin=vmin, vmax=vmax, fontsize=fontsize, postfix=postfix, markersize=markersize,
            figsize=figsize, topowidth=topowidth, preferEright=preferEright, logscale=logscale,
            preferNright=preferNright, saveplot=saveplot, showplot=showplot, dpi=dpi)
            (calcresiduals && doreshist) && plotgausshist(res[i][:], title="Line $(soundings[a].linenum) residuals")
    end
    getphidhist(ϕd; doplot=dophiplot, saveplot=dophiplot, prefix=prefix)
    if calcresiduals
        idxgood = isassigned.(Ref(res), 1:length(res))
        rall = reduce(hcat, res[idxgood])
        doreshist && plotgausshist(vec(rall), title="All Lines residuals")
        map(eachrow(rall)) do r
            n = length(r)
            sqrt(r'r/n) # should return timechannel scaling factors in noise
        end
    end
end

function getphidhist(ϕd; doplot=false, saveplot=false, prefix="", figsize=(8,4), fontsize=11)
    edges=[0,1.1,2, Inf]
    ϕdcounts = fit(Histogram, filter(!isnan, ϕd), edges).weights
    good, bad, ugly = round.(Int, ϕdcounts./sum(ϕdcounts)*100) # % for type of fit in ranges above
    ugly = 100 - (good+bad) # ensure 100%
    if doplot
        fout = prefix == "" ? "summaryfits.png" : prefix*"_summaryfits.png"
        f = figure(figsize=figsize)
        width=[diff(edges[1:end-1]);1]
        bar(edges[1:end-1], [good, bad, ugly], align="edge", width=width)
        xlabel(L"\phi_d")
        ylabel("% of soundings")
        isgood = sum(.!isnan.(ϕd))
        title(prefix*" total soundings: $isgood NaN: $(length(ϕd)-isgood)")
        grid()
        nicenup(f, fsize=fontsize)
        saveplot && savefig(fout, dpi=400)
    end
    good, bad, ugly
end    

function plotconvandlasteachline(soundings, σ, ϕd, delr, delz, resid; 
        idx = nothing, #array of sounding indexes at a line to draw profile
        zall=nothing, cmapσ="turbo", vmin=-2.5, vmax=0.5, fontsize=12,
        figsize=(20,5),
        topowidth=1,
        preferEright = false,
        preferNright = false,
        saveplot = true, 
        showplot = true,
        logscale = false,
        postfix = "",
        yl = nothing,
        markersize = 2,
        dpi = 400)
    @assert !isnothing(zall)
    # ϕd, σ = readingrid(soundings, zall)
    img, gridr, gridz, topofine, R = makegrid(σ, soundings, zall=zall, dz=delz, dr=delr)
    nextra = length(resid) == 0 ? 0 : 1
    height_ratios = nextra == 1 ? [1,1,2,4,0.25] : [1,1,4,0.25]
    fig, ax = plt.subplots(4+nextra, 1, gridspec_kw=Dict("height_ratios" => height_ratios),
        figsize=figsize)
    lname = "Line_$(soundings[1].linenum)_"*postfix
    x0, y0 = soundings[1].X, soundings[1].Y
    if isdefined(soundings[1], :zTx)
        zTx = [s.zTx for s in soundings]
    elseif isdefined(soundings[1], :z_tx)    
        zTx = [s.z_tx for s in soundings] # so we can use for TEMPEST and SPECTREM
    else
        zTx = [s.zTxLM for s in soundings]
    end    
    xend, yend = soundings[end].X, soundings[end].Y
    Z = [s.Z for s in soundings]
    good, bad, ugly = getphidhist(ϕd)
    resminmax = nextra == 0 ? "" : @sprintf(" resid_range: (%.1f,%.1f)", extrema(resid)...)
    fig.suptitle(lname*" Δx=$delr m, Fids: $(length(R)) "*L"\phi_{d_{0-1.1}}:"*" $good "*
        L"\phi_{d_{1.1-2}}:"*" $bad "*L"\phi_{d_{2-\infty}}:"*" $ugly"*resminmax, fontsize=fontsize)
    ax[1].plot(R, ones(length(R)), "--k")
    ax[1].plot(R, ϕd, ".", markersize=markersize)
    ax[1].set_ylim(0.316, maximum(ax[1].get_ylim()))
    ax[1].set_ylabel(L"\phi_d")
    logscale && ax[1].set_yscale("log")
    ax[2].plot(R, zTx)
    ax[2].set_ylabel("zTx m")
    ax[2].sharex(ax[1])
    irow = 2+1
    if nextra == 1 # show resids
        ax[irow].pcolormesh(R, 1:size(resid, 1), resid, vmax=5, vmin=-5, cmap="RdBu_r")
        ax[irow].set_ylabel("window #")
        ax[irow].invert_yaxis()
        ax[irow].sharex(ax[irow-1])
        irow += 1
    end
    imlast = ax[irow].imshow(img, extent=[gridr[1], gridr[end], gridz[end], gridz[1]], cmap=cmapσ, aspect="auto", vmin=vmin, vmax=vmax)
    ax[irow].plot(gridr, topofine, linewidth=topowidth, "-k")
    [a.tick_params(labelbottom=false) for a in ax[1:irow-1]]
    isnothing(idx) || plotprofile(ax[irow], idx, Z, R)
    # eg = extrema(gridr)
    isa(yl, Nothing) || ax[irow].set_ylim(yl...)
    ax[irow].set_ylabel("mAHD")
    ax[irow].set_xlabel("Distance m")
    ax[irow].sharex(ax[irow-1])
    ax[irow].set_xlim(extrema(gridr))
    plotNEWSlabels(soundings, gridr, gridz, [ax[irow]], x0, y0, xend, yend, 
        preferEright=preferEright, preferNright=preferNright; fontsize)
    cb = fig.colorbar(imlast, cax=ax[end], orientation="horizontal")
    cb.set_label("Log₁₀ S/m", labelpad=0)
    nicenup(fig, fsize=fontsize, h_pad=0)
    label = fig._suptitle.get_text()
    VE = round(Int, getVE(ax[end-1]))
    fig.suptitle(label*", VE=$(VE)X"; fontsize)
    saveplot && savefig(lname*".png", dpi=dpi)
    showplot || close(fig)
end    

# driver for gradient based AEM inversion
function loopacrossAEMsoundings(soundings::Array{S, 1}, aem_in, σstart, σ0; 
                            nsequentialiters   =-1,
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
                            minimprovfrac      = nothing,
                            verbose            = false,
                            minimprovkickinstep = round(Int, nstepsmax/2),
                            batchstr           = randstring(15),
                            ) where S<:Sounding

    @assert nsequentialiters  != -1
    nparallelsoundings = nworkers()
    nsoundings = length(soundings)
    zall = zboundarytocenter(aem_in.z[aem_in.nfixed+1:end]) # needed for sounding compression in write
    fname = batchstr*"_"*zipsaveprefix*"_global.log"
    writetogloballog("will require $nsequentialiters iterations of $nparallelsoundings soundings"; iomode="w", fname)
    writetogloballog("starting sequential parallel iterations at $(Dates.now())"; fname)
    for iter = 1:nsequentialiters
        ss = getss_deterministic(iter, nsequentialiters, nparallelsoundings, nsoundings)
        writetogloballog("soundings in loop $iter of $nsequentialiters $ss"; fname)
        pids = workers()
        t2 = time()
        @sync for (i, s) in enumerate(ss)
            aem = makeoperator(aem_in, soundings[s])
            fname_ = soundings[s].sounding_string*"_gradientinv.dat"
            σstart_, σ0_ = map(x->x*ones(length(aem.ρ)-1), [σstart, σ0])
            @async remotecall_fetch(gradientinv, pids[i], σstart_, σ0_, aem;
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
                                                fname              = fname_           ,
                                                minimprovfrac, verbose,
                                                minimprovkickinstep) 
                

        end # @sync
        isfirstparalleliteration = iter == 1 ? true : false
        compresssoundings && compress(soundings[ss[1]:ss[end]], zall, 
            isfirstparalleliteration = isfirstparalleliteration, prefix=zipsaveprefix)
        dt = time() - t2 # seconds
        writetogloballog("done $iter out of $nsequentialiters at $(Dates.now()) in $dt sec"; fname)
    end
end

end # module