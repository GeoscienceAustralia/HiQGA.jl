## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 10001, 8, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using transD_GP
## run McMC
@time transD_GP.main(opt, surface, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
