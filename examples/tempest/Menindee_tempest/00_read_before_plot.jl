## MPI Init same was as on gadi
using MPIClusterManagers, Distributed
import MPI
usempi = false
if usempi
	MPI.Init()
	rank = MPI.Comm_rank(MPI.COMM_WORLD)
	size = MPI.Comm_size(MPI.COMM_WORLD)
	if rank == 0
	    @info "size is $size"
	end
	manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)
	@info "there are $(nworkers()) workers"
	@everywhere @info gethostname()
end
## read the soundings
using PyPlot, Revise, transD_GP
sounding = transD_GP.TEMPEST1DInversion.read_survey_files(
	fname_dat= "Line_2000001",
	fname_specs_halt="electronics_halt.jl",
	frame_height = 19,
	frame_dz = 31,
	frame_dx = 29,
	frame_dy = 30,
	Hxs = [36, 50],
	Hzs = [75, 89],
	Hxp = 69,
	Hzp = 108,
	units = 1e-15,
	yaw_rx = 34,
	pitch_rx = 32,
	roll_rx = 33,
	yaw_tx = 22,
	pitch_tx = 20,
	roll_tx = 21,
	X = 12,
	Y = 13,
	Z = 14,
	fid = 3,
	startfrom = 1,
	skipevery = 10,
	linenum = 1,
	makesounding=true)
## MPI checks
# split into sequential iterations of parallel soundings
nsoundings = length(sounding)
nchainspersounding = 7
ppn = 48
if usempi
	ncores = nworkers()
	@assert mod(ncores+1,nchainspersounding+1) == 0
	@assert mod(ppn, nchainspersounding+1) == 0
	nparallelsoundings = Int((ncores+1)/(nchainspersounding+1))
	nsequentialiters = ceil(Int, nsoundings/nparallelsoundings)
	@info "will require $nsequentialiters iterations of $nparallelsoundings"
end
nsamples = 1001
