## delete this on NCI
# using Distributed
# addprocs(2)
## end delete
## MPI Init
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
transD_GP.loopacrossAEMsoundings(soundings, σstart, σ0,
                            nsequentialiters   = nsequentialiters,              
                            zfixed             = zfixed          ,              
                            ρfixed             = ρfixed          ,              
                            zstart             = zstart          ,              
                            extendfrac         = extendfrac      ,              
                            dz                 = dz              ,              
                            ρbg                = ρbg             ,              
                            nlayers            = nlayers         ,              
                            ntimesperdecade    = ntimesperdecade ,              
                            nfreqsperdecade    = nfreqsperdecade ,              
                            modelprimary       = modelprimary    ,              
                            regtype            = regtype         ,              
                            nstepsmax          = nstepsmax       ,              
                            ntries             = ntries          ,              
                            target             = target          ,              
                            lo                 = lo              ,              
                            hi                 = hi              ,              
                            λ²min              = λ²min           ,              
                            λ²max              = λ²max           ,              
                            λ²frac             = λ²frac          ,              
                            β²                 = β²              ,
                            ntestdivsλ²        = ntestdivsλ²     ,              
                            αmin               = αmin            ,              
                            αmax               = αmax            ,              
                            αfrac              = αfrac           ,              
                            ntestdivsα         = ntestdivsα      ,              
                            regularizeupdate   = regularizeupdate,              
                            knownvalue         = knownvalue      ,              
                            firstvalue         = firstvalue      ,              
                            κ                  = κ               ,              
                            breakonknown       = breakonknown    ,              
                            dobo               = dobo,
                            zipsaveprefix      = zipsaveprefix)              

MPIClusterManagers.stop_main_loop(manager)
rmprocs(workers())
exit()
