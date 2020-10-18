srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
## MPI Init same was as on gadi
using MPIClusterManagers, Distributed
import MPI
usempi = true
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
## info to read data
# datafile
fname_dat = "em_north_line.dat"
# electronics file
fname_specs_halt = "electronics_halt.jl"
# column numbers from hdr file
X, Y = 12, 13
fid = 1
linenum = 2
frame_height = 9
frame_dx = 21
frame_dy = 22
frame_dz = 23
LM_Z = [26,42]
HM_Z = [43,67]
##
using Revise, SkyTEM1DInversion
SkyTEM1DInversion.read_survey_files(fname_dat         = fname_dat,
						             fname_specs_halt = fname_specs_halt,
						             LM_Z             = LM_Z,
									 HM_Z             = HM_Z,
									 frame_height     = frame_height,
									 frame_dz         = frame_dz,
									 frame_dy         = frame_dy,
									 frame_dx         = frame_dx,
									 X                = X,
									 Y                = Y,
									 fid              = fid,
									 linenum          = linenum)
##
sounding = SkyTEM1DInversion.read_survey_files(fname_dat         = fname_dat,
						             fname_specs_halt = fname_specs_halt,
						             LM_Z             = LM_Z,
									 HM_Z             = HM_Z,
									 frame_height     = frame_height,
									 frame_dz         = frame_dz,
									 frame_dy         = frame_dy,
									 frame_dx         = frame_dx,
									 X                = X,
									 Y                = Y,
									 fid              = fid,
									 linenum          = linenum,
									 skipevery        = 10,
									 dotillsounding   = nothing,
									 makesounding     = true)
## MPI checks
# split into sequential iterations of parallel soundings
nsoundings = length(sounding)
nchainspersounding = 5
ppn = 48
if usempi
	ncores = nworkers()
	@assert mod(ncores+1,nchainspersounding+1) == 0
	@assert mod(ppn, nchainspersounding+1) == 0 
	nparallelsoundings = Int((ncores+1)/(nchainspersounding+1))
	nsequentialiters = ceil(Int, nsoundings/nparallelsoundings)
	@info "will require $nsequentialiters iterations of $nparallelsoundings"
end
nsamples = 100001
