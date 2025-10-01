# set up McMC using 4 chains
using Distributed
nsamples, nchains, nchainsatone = 20001, 4, 1
Tmax = 2.50
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
## plot
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
opt.xall[:] = zall
transD_GP.plot_posterior(swd, opt,  burninfrac=0.5, figsize=(4,4), fsize=8, nbins=50, usekde=true, cmappdf = "bone")