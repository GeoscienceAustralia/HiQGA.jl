# MPI Init
using MPIClusterManagers, Distributed
import MPI
MPI.Init()
rank = MPI.Comm_rank(MPI.COMM_WORLD)
size = MPI.Comm_size(MPI.COMM_WORLD)
if rank == 0
    @info "size is $size"
end
manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)
@info "there are $(nworkers()) workers"
@everywhere @info gethostname()
include("01_read_data.jl")
include("02_set_options.jl")
## MPI checks
# split into sequential iterations of parallel soundings
nsoundings = length(soundings)
ncores = nworkers()
nsequentialiters = ceil(Int, nsoundings/ncores)
@info "will require $nsequentialiters iterations of $ncores soundings"
## set up transD_GP
@everywhere using Distributed
@everywhere using HiQGA.transD_GP
## do the parallel soundings
@info "starting"
transD_GP.loopacrossAEMsoundings(soundings, aem, σstart, σ0;
                            nsequentialiters, regtype, nstepsmax, ntries,           
                            lo, hi, λ²min, λ²max, β², knownvalue, breakonknown,              
                            zipsaveprefix)                  

MPIClusterManagers.stop_main_loop(manager)
rmprocs(workers())
exit()
