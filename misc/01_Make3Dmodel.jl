srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Statistics, Distributed, Revise, Tools3D, TransD_GP, MCMC_Driver
##
n, dz, extendfrac  = 50, 1.01, 1.1
z, znall = Tools3D.makezn(n=n, dz=dz, extendfrac=extendfrac)

opt, optdummy, x, y = Tools3D.makeopt(znall=znall, dxfrac = 0.01, dyfrac=0.01)
m = TransD_GP.init(opt)
mdummy = TransD_GP.init(optdummy, m)
Tools3D.slicemodel(m, opt, slicesx=[], slicesy=[], slicesz=[16,19], dz=dz, extendfrac=extendfrac)

ρ = Tools3D.makecubemodel(opt, dz=dz, extendfrac=extendfrac, ρanom=0.2)
slicesx, slicesy, slicesz = [20], [20], [17]
Tools3D.slicemodel(ρ,0,[0. 0.],[0.], opt, slicesx=slicesx, slicesy=slicesy, slicesz=slicesz, dz=dz, extendfrac=extendfrac)
savefig("true_model.png",dpi=300)

d, sd = Tools3D.get_training_data(ρ, opt, dz=dz, extendfrac=extendfrac, zbreak=2.0, fractrain=0.05)
savefig("random_data.png")
Tools3D.calc_simple_RMS(d, ρ, sd)
