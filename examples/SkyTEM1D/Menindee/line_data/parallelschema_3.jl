srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
## load the sounding names
fdataname = ["Line_115651_dbdt_gates_rangeidx_40.mat",
              "Line_115651_dbdt_gates_rangeidx_123.mat",
              "Line_115651_dbdt_gates_rangeidx_205.mat",
              "Line_115651_dbdt_gates_rangeidx_288.mat"]
## same for all
zfixed   = [-1e5]
ρfixed   = [1e12]
zstart = 0.0
extendfrac, dz = 1.06, 2.
nlayers = 40
ρbg = 10.
ntimesperdecade = 10
nfreqsperdecade = 5
## make transD options
using GP
nmin, nmax = 2, 40
K = GP.Mat32()
demean = true
fbounds = [-0.5 2.5]
sdpos = 0.05
sdprop = 0.05
λ, δ = [2], 0.1
save_freq = 25
## split into sequential iterations of parallel soundings
nsoundings = length(fdataname)
ncores = 19
nchainspersounding = 4
@assert mod(ncores+1,nchainspersounding+1) == 0
nparallelsoundings = Int((ncores+1)/(nchainspersounding+1))
nsequentialiters = ceil(Int, nsoundings/nparallelsoundings)
## set up McMC
using Distributed, MCMC_Driver
nsamples, nchains, nchainsatone = 100001, 4, 1
usempi = true
if usempi
    using MPIClusterManagers
    manager = MPIManager(np=ncores)
    addprocs(manager)
else
    addprocs(ncores)
end
@info "there are $(nworkers()) workers"
@everywhere @info gethostname()
Tmax = 2.50
@info "workers are $(workers())"
@everywhere any($srcdir .== LOAD_PATH) || push!(LOAD_PATH, $srcdir)
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
@everywhere using AEM_VMD_HMD, SkyTEM1DInversion
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
            SkyTEM1DInversion.makeoperator(fdataname[s],
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
                                fileprefix = fdataname[s],
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
        #sleep(300)
    end
end
exit()
