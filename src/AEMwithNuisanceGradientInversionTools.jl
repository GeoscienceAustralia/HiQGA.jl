module AEMwithNuisanceGradientInversionTools
using Distributed, Dates, Printf, PyPlot, DelimitedFiles, StatsBase
using ..AbstractOperator, ..CommonToAll, ..GP
import ..AbstractOperator.makeoperator
import ..AbstractOperator.setnuboundsandstartforgradinv
import ..AbstractOperator.Sounding
import ..AbstractOperator.returnforwrite
import ..AbstractOperator.loopacrossAEMsoundings
import ..AbstractOperator.plotconvandlast
import ..gradientinv
export plotconvandlast, loopacrossAEMsoundings
# for deterministic inversions, read in
function compress(soundings, zall, nnu; prefix="", rmfile=true, isfirstparalleliteration=false)
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
        σgrid = vec(A[end,3:end-nnu])
        nu = vec(A[end,end-nnu+1:end])
        for el in [returnforwrite(s)...; vec(zall); σgrid; nu; ϕd]
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
function plotconvandlast(soundings, delr, delz, nufieldnames::Vector{Symbol}; 
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
        lnames = nothing,
        idx = nothing, # array of arrrays per line
        yl = nothing,
        dpi=400)
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)
    nnu = length(nufieldnames)
    nlayers = length(zall)
    fnamecheck = soundings[1].sounding_string*"_gradientinv.dat"
    isfile(fnamecheck) && compress(soundings, zall, nnu) # write everything in one file if not done yet
    fzipped = prefix == "" ? "zipped.dat" : prefix*"_zipped.dat"
    A = readdlm(fzipped)
    σ = A[:,end-nlayers-nnu:end-nnu-1]
    nu = A[:,end-nnu:end-1]
    ϕd = A[:,end]
    for i in 1:nlines
        a = linestartidx[i]
        b = i != nlines ?  linestartidx[i+1]-1 : length(soundings)
        idspec = nothing
        if !isnothing(lnames) # only specific lines wanted, nothing means all lines
            @assert length(lnames) == length(idx)
            doesmatch = findfirst(lnames .== soundings[a].linenum) 
            isnothing(doesmatch) && continue
            @info lnames[doesmatch]
            @show idspec = idx[doesmatch]
        end    
        plotconvandlasteachline(soundings[a:b], view(σ, a:b, :)', view(nu, a:b, :), nufieldnames, view(ϕd, a:b), delr, delz; 
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

function plotconvandlasteachline(soundings, σ, nu, nufieldnames, ϕd, delr, delz; 
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
    nnu = length(nufieldnames)
    img, gridr, gridz, topofine, R = makegrid(σ, soundings, zall=zall, dz=delz, dr=delr)
    fig, ax = plt.subplots(3+nnu, 1, gridspec_kw=Dict("height_ratios" => [1,ones(1+nnu)...,4]),
        figsize=figsize, sharex=true)
    lname = "Line_$(soundings[1].linenum)"*postfix
    x0, y0 = soundings[1].X, soundings[1].Y
    if isdefined(soundings[1], :z_tx) # make this a symbol key TODO
        zTx = [s.z_tx for s in soundings]
    end
    nugiven = reduce(hcat, [map(fldname->[getfield(s, fldname) for s in soundings], nufieldnames)...]) # nsoundings×nnu 
    xend, yend = soundings[end].X, soundings[end].Y
    Z = [s.Z for s in soundings]
    good, bad, ugly = getphidhist(ϕd)
    fig.suptitle(lname*" Δx=$delr m, Fids: $(length(R)) "*L"\phi_{d_{0-1.1}}:"*" $good "*
    L"\phi_{d_{1.1-2}}:"*" $bad "*L"\phi_{d_{2-\infty}}:"*" $ugly", fontsize=fontsize)
    ax = fig.axes
    ax[1].plot(R, ones(length(R)), "--k")
    ax[1].plot(R, ϕd, ".", markersize=markersize)
    ax[1].set_ylim(0.316, maximum(ax[1].get_ylim()))
    ax[1].set_ylabel(L"\phi_d")
    logscale && ax[1].set_yscale("log")
    ax[2].plot(R, zTx)
    ax[2].set_ylabel("zTx m")
    irow = 3
    for (inu, nfname) in enumerate(nufieldnames)
        ax[irow].plot(R, nugiven[:,inu], "--k")
        ax[irow].plot(R, nu[:,inu])
        ax[irow].set_ylabel("$(nfname)")
        irow += 1
    end    
    [a.tick_params(labelbottom=false) for a in ax[1:irow-1]]
    imlast = ax[irow].imshow(img, extent=[gridr[1], gridr[end], gridz[end], gridz[1]], cmap=cmapσ, aspect="auto", vmin=vmin, vmax=vmax)
    ax[irow].plot(gridr, topofine, linewidth=topowidth, "-k")
    isnothing(idx) || plotprofile(ax[irow], idx, Z, R)
    # eg = extrema(gridr)
    isa(yl, Nothing) || ax[3].set_ylim(yl...)
    ax[irow].set_ylabel("mAHD")
    ax[irow].set_xlabel("Distance m")
    fig.colorbar(imlast, ax=ax[irow], location="bottom", shrink=0.6, label="Log₁₀ S/m")
    # ax[irow].set_xlim(extrema(gridr))
    plotNEWSlabels(soundings, gridr, gridz, [ax[irow]], x0, y0, xend, yend, 
    preferEright=preferEright, preferNright=preferNright)
    nicenup(fig, fsize=fontsize)
    label = fig._suptitle.get_text()
    VE = round(Int, getVE(ax[end-1]))
    fig.suptitle(label*", VE=$(VE)X")
    saveplot && savefig(lname*".png", dpi=dpi)
    showplot || close(fig)
end    

# driver for gradient based AEM inversion
function loopacrossAEMsoundings(soundings::Array{S, 1}, aem_in, σstart, σ0, nuboundsΔ; 
                            nsequentialiters   =-1,
                            zstart             = 0.0,
                            extendfrac         = 1.06,
                            dz                 = 2.,
                            nlayers            = 50,
                            regtype            = :R1,
                            nstepsmax          = 10,
                            ntries             = 6,
                            target             = nothing,
                            lo                 = -3.,
                            hi                 = 1.,
                            λ²min              = 0,
                            λ²max              = 8,
                            β²                 = 0.,
                            regularizeupdate   = false,
                            knownvalue         = 0.7,
                            compresssoundings  = true,
                            zipsaveprefix      = "",
                            # optim stuff
                            ntriesnu           = 5,
                            boxiters           = 3, 
                            usebox             = true,
                            reducenuto         = 0.2,
                            debuglevel         = 0,
                            breaknuonknown     = false,       
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
            aem = makeoperator(aem_in, soundings[s])
            nubounds, nustart = setnuboundsandstartforgradinv(aem, nuboundsΔ)
            fname = soundings[s].sounding_string*"_gradientinv.dat"
            σstart_, σ0_ = map(x->x*ones(length(aem.ρ)-1), [σstart, σ0])
            @async remotecall_wait(gradientinv, pids[i], σstart_, σ0_, nustart, aem;
                                                regtype            = regtype         ,              
                                                nstepsmax          = nstepsmax       ,              
                                                ntries             = ntries          ,              
                                                target             = target          ,              
                                                lo                 = lo              ,              
                                                hi                 = hi              ,              
                                                λ²min              = λ²min           ,              
                                                λ²max              = λ²max           ,              
                                                β²                 = β²              ,
                                                regularizeupdate   = regularizeupdate,              
                                                knownvalue         = knownvalue      ,              
                                                fname              = fname           ,
                                                ntriesnu, nubounds, boxiters, usebox, reducenuto, debuglevel, breaknuonknown) 
                

        end # @sync
        isfirstparalleliteration = iter == 1 ? true : false
        compresssoundings && compress(soundings[ss[1]:ss[end]], zall, size(nuboundsΔ, 1),
            isfirstparalleliteration = isfirstparalleliteration, prefix=zipsaveprefix)
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end

end # module