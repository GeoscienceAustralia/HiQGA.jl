## make options for the stationary GP
nmin, nmax = 2, 40
pnorm = 2.
K = transD_GP.GP.OrstUhn()
demean = false
fbounds = [-0.5 3.25]
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.1diff(fbounds, dims=2)[:]
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
λ, δ = [2], 0.1
fdataname = "nutest_"
sampledc = true
sdev_dc = 0.012*diff(fbounds, dims=2)[:]
## Initialize a stationary GP using these options
Random.seed!(12)
opt = transD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        fdataname = fdataname,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        sampledc = sampledc,
                        sdev_dc = sdev_dc,
                        pnorm = pnorm,
                        save_freq = 25,
                        quasimultid = false,
                        K = K,
                        updatenuisances = true
                        )
## nuisance options
# we don't treat *every* geometry param as nuisance.
# this is achieved by having the fixed parameters set as nuisances with
# zero-width bounds. These values are *fractions*
optn = transD_GP.OptionsNuisance(opt;
    sdev = [
        0.0 #zTx
        0.005 #zRx
        0.005 #x_rx
        0.0 #y_rx
        0.0 #rx_roll
        0.0 #rx_pitch
        0.0 #rx_yaw
        0.0 #tx_roll
        0.0 #tx_pitch
        0.0 #tx_yaw
           ],
    bounds = [
      zTx               zTx
      zRx - 3.0         zRx + 3.0
      x_rx - 3.0        x_rx + 3.0
      y_rx              y_rx
      rx_roll           rx_roll
      rx_pitch          rx_pitch
      rx_yaw            rx_yaw
      tx_roll           tx_roll
      tx_pitch          tx_pitch
      tx_yaw            tx_yaw
          ],
    updatenuisances = true,
        )
