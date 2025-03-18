## MPI Init
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
include("02_set_options.jl")
## set up McMC
@everywhere using Distributed
@everywhere using HiQGA.transD_GP
## do the parallel soundings
dr = 15
yl = []
idx = []
qp1, qp2 = 0.1, 0.9
burninfrac = 0.25
vmin, vmax = extrema(-fbounds)
using Dates, Printf
dt = time()
@info "starting at $(Dates.now())"
transD_GP.summaryAEMimages(soundings, opt; qp1, qp2,
    zall, dr, dz, cmap="turbo", burninfrac, figsize=(16,12), yl,
    vmin=vmin, vmax=vmax, idx, useML, fontsize=14,
    preferEright=true, saveplot=true, showplot = false,)

@info "done! at $(Dates.now()) in "*@sprintf("%.2f seconds", time() - dt)
## stop and exit
MPIClusterManagers.stop_main_loop(manager)
rmprocs(workers())
exit()
