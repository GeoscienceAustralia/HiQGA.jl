## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 2001, 7, 1
Tmax = 2.5
addprocs(nchains)
##init packages on workers
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using HiQGA.transD_GP
## run McMC
@time transD_GP.main(opt, optn, tempest, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)

## clean up
rmprocs(workers())
