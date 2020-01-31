nsamples, nchains, nchainsatone = 200001, 16, 1
Tmax = 2.5 

addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(opt, d, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)

