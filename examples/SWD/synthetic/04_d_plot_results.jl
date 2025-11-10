## plot
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
opt.xall[:] = zall
transD_GP.plot_posterior(swd, opt,  burninfrac=0.25, figsize=(4,6), 
    fsize=8, nbins=50, usekde=true, cmappdf = "bone_r",
    qp1=0.1, qp2=0.9)
ax = gcf().axes
ax[1].step(vs, zboundaries, "b--")