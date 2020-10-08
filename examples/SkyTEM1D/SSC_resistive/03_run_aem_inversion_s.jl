##MPI stuff
using MPIClusterManagers, Distributed
import MPI 

MPI.Init()
rank = MPI.Comm_rank(MPI.COMM_WORLD)
size = MPI.Comm_size(MPI.COMM_WORLD)

manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)
@everywhere @info gethostname()

## set up McMC
nsamples, nchains, nchainsatone = 100001, 4, 1
Tmax = 2.50
@info "workers are $(workers())"
@everywhere any($srcdir .== LOAD_PATH) || push!(LOAD_PATH, $srcdir)
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(opt, optdummy, aem, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
MPIClusterManagers.stop_main_loop(manager)
exit()

