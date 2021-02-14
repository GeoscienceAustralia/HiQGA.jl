using transD_GP, Distributed
## make options for the stationary GP
nmin, nmax = 2, 40
pnorm = 2.
K = transD_GP.GP.Mat32()
demean = true
fbounds = [-0.5 2.5]
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.05diff(fbounds, dims=2)[:]
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
updatenonstat = false
needλ²fromlog = false
λ, δ = [2], 0.1

## Initialize a stationary GP using these options
Random.seed!(12)
opt = transD_GP.OptionsStat(nmin = nmin,
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
                        save_freq = 25,
                        quasimultid = false,
                        K = K,
                        needλ²fromlog = needλ²fromlog,
                        updatenonstat = updatenonstat,
                        updatenuisances = true
                        )
## Initialize options for the dummy nonstationary properties GP
Random.seed!(13)
optdummy = transD_GP.OptionsNonstat(opt,
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

## nuisance options
#we don't treat *every* geometry param as nuisance.
#this is achieved with some computational waste by
#having the fixed parameters set as nuisances with
#zero-width bounds.
optn = transD_GP.OptionsNuisance(opt)

optn.sdev =
    [
        0.0 #zTx
        1.0 #zRx
        1.0 #x_rx
        0.0 #y_rx
        0.0 #rx_roll
        0.4 #rx_pitch
        0.0 #rx_yaw
        0.0 #tx_roll
        0.0 #tx_pitch
        0.0 #tx_yaw
    ]
optn.bounds =
    [
      zTx               zTx
      zRx - 5.0         zRx + 5.0
      x_rx - 5.0        x_rx + 5.0
      y_rx              y_rx
      rx_roll           rx_roll
      rx_pitch - 2.0    rx_pitch + 2.0
      rx_yaw            rx_yaw
      tx_roll           tx_roll
      tx_pitch          tx_pitch
      tx_yaw            tx_yaw
    ]
optn.nnu = 10

optn.updatenuisances = true
##
optn


##debug stuff
mn = transD_GP.init(optn)
mstat = transD_GP.init(opt)
mnstat = transD_GP.init(optdummy, mstat)

##
initmisfit = transD_GP.get_misfit(mstat, opt, tempest)
altinit = transD_GP.get_misfit(mstat, mn, opt, tempest)

statn = transD_GP.Stats()

fpc = open("debugcosts.bin", "w")
fpv = open("debugvals.bin", "w")
wpn = transD_GP.Writepointers_nuisance(fpc,fpv)

Juno.@enter transD_GP.do_mcmc_step(mn, mstat, mnstat, optn, statn, [initmisfit],
    tempest, 1.0, 1, wpn)
