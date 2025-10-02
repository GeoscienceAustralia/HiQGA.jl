# set up McMC using 4 chains
using Distributed
nsamples, nchains, nchainsatone = 50001, 6, 1
Tmax = 2.2
addprocs(nchains)
@info "workers are $(workers())"
@everywhere begin
    using Distributed
    using HiQGA.transD_GP
    include("SWDPhysics.jl")
    include("SWDInversion.jl")
    using .SWDInversion
end    
## run McMC - sample a uniform prior between 0 and 1 at depths between 0 and 1 with 5 nuclei
@time transD_GP.main(opt, swd, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())