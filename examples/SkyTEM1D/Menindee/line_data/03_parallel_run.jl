srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
## set up McMC
@everywhere any($srcdir .== LOAD_PATH) || push!(LOAD_PATH, $srcdir)
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
@everywhere using SkyTEM1DInversion
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
            SkyTEM1DInversion.makeoperator(sounding[s],
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
        r_opt_and_optdummy = @spawnat pids[2] begin
            SkyTEM1DInversion.make_tdgp_statmode_opt(znall = fetch(r_aem_and_znall)[2],
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
                    MCMC_Driver.main(fetch(r_opt_and_optdummy)...,
                    fetch(r_aem_and_znall)[1],
                    collect(pids[2:end]),
                    Tmax=Tmax,
                    nsamples=nsamples,
                    nchainsatone=nchainsatone)
        end
        @info "done $i"
    end
end
exit()
