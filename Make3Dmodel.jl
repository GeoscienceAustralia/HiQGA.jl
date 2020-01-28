any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using PyPlot, Distributed, Revise, Tools3D, TransD_GP, MCMC_Driver
##
n, dz, extendfrac  = 20, 1.6, 1.245
z, znall = Tools3D.makezn(n=n, dz=dz, extendfrac=extendfrac)

opt, x, y = Tools3D.makeopt(znall=znall)
m = TransD_GP.init(opt)
Tools3D.slicemodel(m, opt, slicesx=[], slicesy=[], slicesz=[16,19], dz=dz, extendfrac=extendfrac)

ρ = Tools3D.makecubemodel(opt, dz=dz, extendfrac=extendfrac, ρanom=0.2)
Tools3D.slicemodel(ρ,0,[0. 0.],[0.], opt, slicesx=[20], slicesy=[], slicesz=[], dz=dz, extendfrac=extendfrac)
Tools3D.slicemodel(ρ,0,[0. 0.],[0.], opt, slicesx=[], slicesy=[20], slicesz=[], dz=dz, extendfrac=extendfrac)
Tools3D.slicemodel(ρ,0,[0. 0.],[0.], opt, slicesx=[], slicesy=[], slicesz=[17], dz=dz, extendfrac=extendfrac)

d, sd = Tools3D.get_training_data(ρ, opt, dz=dz, extendfrac=extendfrac, zbreak=100.0, fractrain=0.05)

Tools3D.calc_simple_RMS(d, ρ, sd)
## set up McMC
nsamples, nchains, nchainsatone = 500001, 8, 1
Tmax = 2.5

addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(opt, d, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
##
MCMC_Driver.getchi2forall(opt)
Tools3D.plot_last_target_model(opt, d, sd; slicesx=[20], slicesy=[], slicesz=[], dz=dz, extendfrac=extendfrac)
Tools3D.plot_last_target_model(opt, d, sd; slicesx=[], slicesy=[20], slicesz=[], dz=dz, extendfrac=extendfrac)
Tools3D.plot_last_target_model(opt, d, sd; slicesx=[], slicesy=[], slicesz=[17], dz=dz, extendfrac=extendfrac)
rmprocs(workers())
