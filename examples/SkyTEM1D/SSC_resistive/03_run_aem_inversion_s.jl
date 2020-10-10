## MPI stuff and parallel tempering chains
using MPIClusterManagers, Distributed
nchains = 4
manager = MPIManager(np=nchains)
addprocs(manager)
@info "there are $(nworkers()) workers"
@everywhere @info gethostname()
## set up McMC
nsamples, nchainsatone = 100001, 1
Tmax = 2.50
@info "workers are $(workers())"
@everywhere any($srcdir .== LOAD_PATH) || push!(LOAD_PATH, $srcdir)
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(opt, optdummy, aem, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
exit()
