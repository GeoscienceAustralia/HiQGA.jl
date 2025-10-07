# set up McMC using 4 chains
using Distributed
nsamples, nchains, nchainsatone = 50001, 4, 1 # make sure your computer has nchains+1 cpus!
Tmax = 2.5
addprocs(nchains)
@info "workers are $(workers())"
@everywhere begin
    using Distributed
    using HiQGA.transD_GP
    include("../src/SWDPhysics.jl")
    include("../src/SWDInversion.jl")
    using .SWDInversion
end    
## run McMC - sample a uniform prior between 0 and 1 at depths between 0 and 1 with 5 nuclei
opt.dispstatstoscreen = true # if you don't want acceptance stats on your terminal set false
@time transD_GP.main(opt, swd, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())