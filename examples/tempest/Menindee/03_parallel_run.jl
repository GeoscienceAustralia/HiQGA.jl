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
@assert mod(ppn, nchainspersounding+1) == 0
nparallelsoundings = Int((ncores+1)/(nchainspersounding+1))
nsequentialiters = ceil(Int, nsoundings/nparallelsoundings)
@info "will require $nsequentialiters iterations of $nparallelsoundings"
## set up McMC
@everywhere using Distributed
@everywhere using transD_GP
## do the parallel soundings
@info "starting"
transD_GP.TEMPEST1DInversion.loopacrosssoundings(soundings;
                    nsequentialiters   = nsequentialiters,
                    nparallelsoundings = nparallelsoundings,
                    zfixed             = zfixed,
                    ρfixed             = ρfixed,
                    zstart             = zstart,
                    extendfrac         = extendfrac,
                    dz                 = dz,
                    ρbg                = ρbg,
                    nlayers            = nlayers,
                    ntimesperdecade    = ntimesperdecade,
                    nfreqsperdecade    = nfreqsperdecade,
                    Tmax               = Tmax,
                    nsamples           = nsamples,
                    nchainsatone       = nchainsatone,
                    nchainspersounding = nchainspersounding,
                    nmin               = nmin,
                    nmax               = nmax,
                    K                  = K,
                    demean             = demean,
                    sampledc           = sampledc,
                    sddc               = sddc,
                    sdpos              = sdpos,
                    sdprop             = sdprop,
                    fbounds            = fbounds,
                    save_freq          = save_freq,
                    λ                  = λ,
                    δ                  = δ,
                    useML              = useML,   
                    nuisance_bounds    = nuisance_bounds,
                    nuisance_sdev      = nuisance_sdev,
                    updatenuisances    = updatenuisances,
                    vectorsum          = vectorsum,
                    dispstatstoscreen  = false)


MPIClusterManagers.stop_main_loop(manager)
rmprocs(workers())
exit()
