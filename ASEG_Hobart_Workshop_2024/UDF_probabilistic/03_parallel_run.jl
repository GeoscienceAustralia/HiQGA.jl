using MPIClusterManagers, Distributed, Dates
import MPI
MPI.Init()
rank = MPI.Comm_rank(MPI.COMM_WORLD)
sz = MPI.Comm_size(MPI.COMM_WORLD)
if rank == 0
    @info "size is $sz"
end
manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)
@info "there are $(nworkers()) workers at $(Dates.now())"
include("01_read_data.jl")
include("02_set_options.jl")
## set up McMC
@everywhere using Distributed
@everywhere using HiQGA.transD_GP
## do the parallel soundings
@info "starting at $(Dates.now())"
transD_GP.loopacrossAEMsoundings(soundings, aem, opt;
                    Tmax, nsamples, nchainsatone, nchainspersounding, ppn)

MPIClusterManagers.stop_main_loop(manager)
rmprocs(workers())
exit()
