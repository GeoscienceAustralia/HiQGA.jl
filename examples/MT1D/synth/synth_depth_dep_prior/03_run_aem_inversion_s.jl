## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 300001, 4, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using HiQGA.transD_GP
## run McMC
@time transD_GP.main(opt, F, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
## close the worker pool
rmprocs(workers())
