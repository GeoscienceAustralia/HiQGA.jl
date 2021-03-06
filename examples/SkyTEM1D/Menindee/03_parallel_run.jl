## set up McMC
@everywhere using Distributed
@everywhere using transD_GP
## do the parallel soundings
@info "starting"
for iter = 1:nsequentialiters
    if iter<nsequentialiters
        @show ss = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
    else
        @show ss = (iter-1)*nparallelsoundings+1:nsoundings
    end
    @sync for (i, s) in enumerate(ss)
        @show pids = (i-1)*nchainspersounding+i:i*(nchainspersounding+1)
        r_aem_and_znall = @spawnat pids[2] begin
            transD_GP.SkyTEM1DInversion.makeoperator(sounding[s],
                               zfixed = zfixed,
                               ρfixed = ρfixed,
                               zstart = zstart,
                               extendfrac = extendfrac,
                               dz = dz,
                               ρbg = ρbg,
                               nlayers = nlayers,
                               ntimesperdecade = ntimesperdecade,
                               nfreqsperdecade = nfreqsperdecade,
                               showgeomplot = false,
                               plotfield = false)
        end
        r_opt = @spawnat pids[2] begin
            transD_GP.SkyTEM1DInversion.make_tdgp_opt(znall = fetch(r_aem_and_znall)[2],
                                fileprefix = sounding[s].sounding_string,
                                nmin = nmin,
                                nmax = nmax,
                                K = K,
                                demean = demean,
                                sdpos = sdpos,
                                sdprop = sdprop,
                                fbounds = fbounds,
                                save_freq = save_freq,
                                λ = λ,
                                δ = δ,
                                )
        end
        @spawnat pids[1] begin
                    transD_GP.main(fetch(r_opt),
                    fetch(r_aem_and_znall)[1],
                    collect(pids[2:end]),
                    Tmax=Tmax,
                    nsamples=nsamples,
                    nchainsatone=nchainsatone)
        end
    end
    @info "done $iter out of $nsequentialiters"
end
MPIClusterManagers.stop_main_loop(manager)
rmprocs(workers())
exit()
