using Distributed
addprocs(7) # number of parallel cores
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
transD_GP.SkyTEM1DInversion.loopacrosssoundings(soundings, σstart, σ0,
                            nsequentialiters   = nsequentialiters,              
                            zfixed             = zfixed          ,              
                            ρfixed             = ρfixed          ,              
                            zstart             = zstart          ,              
                            extendfrac         = extendfrac      ,              
                            dz                 = dz              ,              
                            ρbg                = ρbg             ,              
                            nlayers            = nlayers         ,                                 
                            regtype            = regtype         ,              
                            nstepsmax          = nstepsmax       ,              
                            ntries             = ntries          ,              
                            target             = target          ,              
                            lo                 = lo              ,              
                            hi                 = hi              ,              
                            λ²min              = λ²min           ,              
                            λ²max              = λ²max           ,              
                            β²                 = β²              ,
                            knownvalue         = knownvalue      ,              
                            breakonknown       = breakonknown    ,              
                            )                  
rmprocs(workers())
