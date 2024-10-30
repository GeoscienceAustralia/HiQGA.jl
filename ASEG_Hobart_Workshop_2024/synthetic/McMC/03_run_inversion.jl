# # Set up McMC
# Specify parallel tempering on multiple workers
nsamples, nchains, nchainsatone = 2001, 5, 1
Tmax = 2.50 # maximum annealing temperature
##
# # Add parallel workers
using Distributed
nprocs() > 1 && rmprocs(workers()) # remove workers from earlier run
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using HiQGA.transD_GP
##
# # Run McMC with options on multiple workers
@time transD_GP.main(opt, aem, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
