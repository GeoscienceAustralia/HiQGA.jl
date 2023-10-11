module AEMnoNuisanceGradientInversionTools
using Distributed, Dates, Printf, PyPlot, DelimitedFiles, StatsBase
using ..AbstractOperator, ..CommonToAll, ..GP
import ..AbstractOperator.makeoperator
import ..AbstractOperator.Sounding
import ..AbstractOperator.returnforwrite
import ..AbstractOperator.loopacrossAEMsoundings
import ..AbstractOperator.plotconvandlast
import ..AbstractOperator.plotmodelfield!
import ..gradientinv
export plotconvandlast, loopacrossAEMsoundings
# for deterministic inversions, read in
function compress(soundings, zall; prefix="", rmfile=true, isfirstparalleliteration=false)
    fname = soundings[1].sounding_string*"_gradientinv.dat"    
    !isfile(fname) && throw(AssertionError("file does not exist perhaps soundings already zipped?"))
    fout = prefix == "" ? "zipped.dat" : prefix*"_zipped.dat"
    iomode = "w"
    if isfile(fout)
        isfirstparalleliteration && throw(AssertionError("Zipped file "*fout*" exists, will not overwrite!"))
        iomode = "a"
    end    
    io = open(fout, iomode)
    for (i, s) in enumerate(soundings)
        fname = s.sounding_string*"_gradientinv.dat"
        A = readdlm(fname)
        ϕd = A[end,2]
        σgrid = vec(A[end,3:end])
        for el in [returnforwrite(s)...; vec(zall); σgrid; ϕd]
            msg = @sprintf("%.4f\t", el)
            write(io, msg)
        end
        write(io, "\n")                
        flush(io) # slower but ensures write is complete
        rmfile && rm(fname)
    end
    close(io)    
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
        dpi=400)
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)
    nlayers = length(zall)
    fnamecheck = soundings[1].sounding_string*"_gradientinv.dat"
    isfile(fnamecheck) && compress(soundings, zall) # write everything in one file if not done yet
    fzipped = prefix == "" ? "zipped.dat" : prefix*"_zipped.dat"
    A = readdlm(fzipped)
    σ = A[:,end-nlayers:end-1]
    ϕd = A[:,end]
    for i in 1:nlines
        a, b = linestartend(linestartidx, i, nlines, soundings)
        continueflag, idspec = docontinue(lnames, idx, soundings, a, b)
        continueflag && continue
        if plotforward && !isempty(idspec) && !isnothing(aem_in)
            for id in idspec
                aem = makeoperator(aem_in, soundings[a:b][id])
                m = -vec(σ[a:b,:][id,:]) #log 10 ρ
                plotmodelfield!(aem, m)
                gcf().suptitle("Line $(soundings[a].linenum) index:$id")
                nicenup(gcf())
            end    
        end   
        plotconvandlasteachline(soundings[a:b], view(σ, a:b, :)', view(ϕd, a:b), delr, delz; 
            zall = zall, idx=idspec, yl=yl,
            cmapσ=cmapσ, vmin=vmin, vmax=vmax, fontsize=fontsize, postfix=postfix, markersize=markersize,
            figsize=figsize, topowidth=topowidth, preferEright=preferEright, logscale=logscale,
            preferNright=preferNright, saveplot=saveplot, showplot=showplot, dpi=dpi)
    end
    getphidhist(ϕd, doplot=true, saveplot=saveplot, prefix=prefix)
    nothing
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

function plotconvandlasteachline(soundings, σ, ϕd, delr, delz; 
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
    fig, ax = plt.subplots(4, 1, gridspec_kw=Dict("height_ratios" => [1,1,4,0.25]),
        figsize=figsize)
    lname = "Line_$(soundings[1].linenum)"*postfix
    x0, y0 = soundings[1].X, soundings[1].Y
    if isdefined(soundings[1], :zTx)
        zTx = [s.zTx for s in soundings]
    else
        zTx = [s.zTxLM for s in soundings]
    end    
    xend, yend = soundings[end].X, soundings[end].Y
    Z = [s.Z for s in soundings]
    good, bad, ugly = getphidhist(ϕd)
    fig.suptitle(lname*" Δx=$delr m, Fids: $(length(R)) "*L"\phi_{d_{0-1.1}}:"*" $good "*
    L"\phi_{d_{1.1-2}}:"*" $bad "*L"\phi_{d_{2-\infty}}:"*" $ugly", fontsize=fontsize)
    ax[1].plot(R, ones(length(R)), "--k")
    ax[1].plot(R, ϕd, ".", markersize=markersize)
    ax[1].set_ylim(0.316, maximum(ax[1].get_ylim()))
    ax[1].set_ylabel(L"\phi_d")
    logscale && ax[1].set_yscale("log")
    ax[2].plot(R, zTx)
    ax[2].set_ylabel("zTx m")
    [a.tick_params(labelbottom=false) for a in ax[1:2]]
    ax[2].sharex(ax[1])
    imlast = ax[3].imshow(img, extent=[gridr[1], gridr[end], gridz[end], gridz[1]], cmap=cmapσ, aspect="auto", vmin=vmin, vmax=vmax)
    ax[3].plot(gridr, topofine, linewidth=topowidth, "-k")
    isnothing(idx) || plotprofile(ax[3], idx, Z, R)
    # eg = extrema(gridr)
    isa(yl, Nothing) || ax[3].set_ylim(yl...)
    ax[3].set_ylabel("mAHD")
    ax[3].set_xlabel("Distance m")
    ax[3].sharex(ax[2])
    ax[3].set_xlim(extrema(gridr))
    plotNEWSlabels(soundings, gridr, gridz, [ax[3]], x0, y0, xend, yend, 
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
                            ) where S<:Sounding

    @assert nsequentialiters  != -1
    nparallelsoundings = nworkers()
    nsoundings = length(soundings)
    zall = zboundarytocenter(aem_in.z[aem_in.nfixed+1:end]; fudgelast=false) # needed for sounding compression in write
    for iter = 1:nsequentialiters
        if iter<nsequentialiters
            ss = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
        else
            ss = (iter-1)*nparallelsoundings+1:nsoundings
        end
        @info "soundings in loop $iter of $nsequentialiters", ss
        pids = workers()
        @sync for (i, s) in enumerate(ss)
            aem = makeoperator(aem_in, soundings[s])
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

end # module