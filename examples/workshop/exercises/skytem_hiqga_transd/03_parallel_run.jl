using Distributed
# paralleization checks
addprocs(nchainspersounding)
nsoundings = length(soundings)
ncores = nworkers()
@assert mod(ncores+1,nchainspersounding+1) == 0
@assert mod(ppn, nchainspersounding+1) == 0
nparallelsoundings = Int((ncores+1)/(nchainspersounding+1))
nsequentialiters = ceil(Int, nsoundings/nparallelsoundings)
@info "will require $nsequentialiters iterations of $nparallelsoundings parallel soundings"
## set up McMC
@everywhere using Distributed
@everywhere using HiQGA.transD_GP
## do the parallel soundings
@info "starting"
transD_GP.SkyTEM1DInversion.loopacrosssoundings(soundings, opt;
                    nsequentialiters   = nsequentialiters,
                    nparallelsoundings = nparallelsoundings,
                    zfixed             = zfixed,
                    ρfixed             = ρfixed,
                    zstart             = zstart,
                    extendfrac         = extendfrac,
                    dz                 = dz,
                    ρbg                = ρbg,
                    nlayers            = nlayers,
                    Tmax               = Tmax,
                    nsamples           = nsamples,
                    nchainspersounding = nchainspersounding)

rmprocs(workers())