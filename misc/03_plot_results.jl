MCMC_Driver.getchi2forall(opt)
savefig("convergence.png", dpi=300)

irange=1:10:101
Tools3D.savepngs(opt, irange, dz=dz, extendfrac=extendfrac, slicesx=slicesx, slicesy=slicesy, slicesz=slicesz)

M  = Tools3D.assembleTat1(opt;burninfrac=0.33)
Mmean = mean(M, dims=2)
Msd = sqrt.(var(M, dims=2))

Tools3D.slicemodel(vec(Mmean),0,[0. 0.],[0.], opt, slicesx=slicesx, slicesy=slicesy, slicesz=slicesz, dz=dz, extendfrac=extendfrac)
savefig("meanmodel.png", dpi=300)
Tools3D.slicemodel(vec(Msd),0,[0. 0.],[0.], opt, slicesx=slicesx, slicesy=slicesy, slicesz=slicesz, dz=dz, extendfrac=extendfrac)
savefig("sdevmodel.png")
