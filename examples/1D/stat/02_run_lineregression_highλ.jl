## make options for the stationary GP
λ = [.1]
δ = 0.05
nmin, nmax = 2, 30
pnorm = 2.
K = GP.Mat32()
demean = true
fbounds = permutedims([extrema(ynoisy[.!isnan.(ynoisy)])...])
sdev_pos = [0.05abs(diff([extrema(x)...])[1])]
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
xall = permutedims(collect(x))
xbounds = permutedims([extrema(x)...])
updatenonstat = false
needλ²fromlog = false
## Initialize a stationary GP using these options
Random.seed!(12)
opt = TransD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K,
                        timesλ = 3.6,
                        needλ²fromlog = needλ²fromlog,
                        updatenonstat = updatenonstat
                        )
## Initialize options for the dummy nonstationary properties GP
Random.seed!(13)
optdummy = TransD_GP.OptionsNonstat(opt,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K
                        )
## set up McMC
nsamples, nchains, nchainsatone = 200001, 4, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere any($srcdir .== LOAD_PATH) || push!(LOAD_PATH, $srcdir)
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(opt, optdummy, line, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
GeophysOperator.getchi2forall(opt, fsize=8)
ax = gcf().axes;
ndata = sum(.!isnan.(line.d[:]))
ax[3].set_ylim(ndata/2 - 20, ndata/2 + 20)
ax[4].set_ylim(ndata/2 - 20, ndata/2 + 20)
savefig("line_conv_s.png", dpi=300)
GeophysOperator.plot_posterior(line, opt,
    burninfrac=0.5, figsize=(4,4), fsize=8, nbins=50)
ax = gcf().axes
ax[1].plot(ynoisy, x, ".w", alpha=0.1, markersize=10)
ax[1].plot(y, x, "--w", alpha=0.5)
ax[1].set_xlim(fbounds...)
savefig("jump1D_high.png", dpi=300)
