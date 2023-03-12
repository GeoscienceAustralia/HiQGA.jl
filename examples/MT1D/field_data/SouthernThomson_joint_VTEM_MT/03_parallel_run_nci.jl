using MPIClusterManagers, Distributed
import MPI
MPI.Init()
rank = MPI.Comm_rank(MPI.COMM_WORLD)
sz = MPI.Comm_size(MPI.COMM_WORLD)
if rank == 0
    @info "size is $sz"
end
manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)
@info "there are $(nworkers()) workers"
include("01_read_data.jl")
include("02_set_aem_options.jl")
## set up McMC
@everywhere using Distributed
@everywhere using HiQGA.transD_GP
if useML
    # aem.useML = true
    FMT.useML = true
    opt.fdataname*="withML_"
end  
nsamples, nchains, nchainsatone = 500_001, nworkers(), 1
Tmax = 2.50
F = transD_GP.GenericJointInversion.JointOperator([aem, FMT])
## run McMC
@time transD_GP.main(opt, F, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
MPIClusterManagers.stop_main_loop(manager)
rmprocs(workers())
exit()
