# set up McMC using 4 chains
using Distributed
nsamples, nchains, nchainsatone = 20001, 4, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere begin
    using Distributed
    using HiQGA.transD_GP
    include("BarPhysics.jl")
    include("BarPhysicsInversion.jl")
    using .BarPhysicsInversion
end    
## run McMC
@time transD_GP.main(opt, FOp, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
transD_GP.plot_posterior(FOp, opt,  burninfrac=0.5, figsize=(8.5,4), fsize=8, nbins=50)