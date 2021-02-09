## set up McMC
nsamples, nchains, nchainsatone = 20001, 4, 1
Tmax = 2.50
addprocs(nchains)
##init packages on workers
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using transD_GP
## run McMC
@time transD_GP.main(opt, optdummy, optn, tempest, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)

## clean up
rmprocs(workers())
