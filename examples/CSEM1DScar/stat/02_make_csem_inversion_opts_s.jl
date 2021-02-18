using GP, TransD_GP, GeophysOperator, MCMC_Driver, Distributed
## make options for the stationary GP
nmin, nmax = 2, 40
pnorm = 2.
K = GP.Mat32()
demean = true
fbounds = [-0.5 2.3]
sdev_pos = [0.025abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.02*diff(fbounds, dims=2)[:]
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
updatenonstat = false
needλ²fromlog = false
λ, δ = [2], 0.1
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
                        timesλ = 3.,
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