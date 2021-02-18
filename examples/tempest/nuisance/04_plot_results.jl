transD_GP.getchi2forall(optn)
gcf()
##
opt.xall[:] .= zall
transD_GP.plot_posterior(tempest, opt, burninfrac=0.5, figsize=(4,4),
    cmappdf="inferno", qp1=0.01, qp2=0.99)
step(log10.(ρ[2:end]), z[2:end], color="w", alpha=0.4)
gcf()

## nuisance histograms
transD_GP.plot_posterior(tempest, optn, burninfrac=0.5, nbins=100)
gcf()
