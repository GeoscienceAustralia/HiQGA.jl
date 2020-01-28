any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using Revise, Tools3D, TransD_GP, PyPlot
##
n, dz, extendfrac  = 20, 1.6, 1.245
z, znall = Tools3D.makezn(n=n, dz=dz, extendfrac=extendfrac)

opt, x, y = Tools3D.makeopt(znall=znall)
m = TransD_GP.init(opt)
Tools3D.slicemodel(m, opt, slicesx=[], slicesy=[], slicesz=[16,19], dz=dz, extendfrac=extendfrac)
Tools3D.slicemodel(m, opt, slicesx=[30], slicesy=[], slicesz=[], dz=dz, extendfrac=extendfrac)

ρ = Tools3D.makecubemodel(opt, dz=dz, extendfrac=extendfrac, ρanom=0.2)
Tools3D.slicemodel(ρ,0,[0. 0.],[0.], opt, slicesx=[20], slicesy=[], slicesz=[], dz=dz, extendfrac=extendfrac)
Tools3D.slicemodel(ρ,0,[0. 0.],[0.], opt, slicesx=[], slicesy=[20], slicesz=[], dz=dz, extendfrac=extendfrac)
Tools3D.slicemodel(ρ,0,[0. 0.],[0.], opt, slicesx=[], slicesy=[], slicesz=[17], dz=dz, extendfrac=extendfrac)

d = Tools3D.get_training_data(ρ, opt, dz=dz, extendfrac=extendfrac, zbreak=100.0, fractrain=0.05)
