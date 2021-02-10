## model
zfixed   = [-1e6,   0, 40, 60]
rho      = [1e12,   20, 1, 100]
# operator
ntimesperdecade = 10
nfreqsperdecade = 5
nkᵣeval = 60
zTx = -120
zRx = -80
x_rx = -115.
y_rx = 0.
rx_roll = 0.
rx_pitch = 0.
rx_yaw = 0.
tx_roll = 10. #ACHTUNG!
tx_pitch = 10. #ACHTUNG!
tx_yaw = 0.
addprimary = true
tempest = transD_GP.TEMPEST1DInversion.Bfield(
  				     ntimesperdecade = ntimesperdecade,
                 nfreqsperdecade = nfreqsperdecade,
					  nkᵣeval = nkᵣeval,
                 times  = times,
                 ramp   = ramp,
					  zTx    = zTx,
					  zRx    = zRx,
					  x_rx   = x_rx,
					  y_rx   = y_rx,
					  rx_roll = rx_roll,
					  rx_pitch = rx_pitch,
					  rx_yaw = rx_yaw,
					  tx_roll = tx_roll,
					  tx_pitch = tx_pitch,
					  tx_yaw = tx_yaw,
					  addprimary = addprimary)
# X and Z components in fT from GA-AEM
ross_primary = [33.0699      -13.689]
ross_secondary = [
6.148033      -8.96926
   4.254406     -7.237281
   3.515803     -6.429036
   3.037861     -5.838994
   2.680435      -5.36655
   2.369648     -4.949131
   2.011034      -4.45433
    1.59146     -3.831225
   1.147514     -3.094577
  0.7279429     -2.289014
  0.3895116      -1.50605
  0.1698442    -0.8633888
0.05873083    -0.4281239
0.01569159    -0.1862741
0.002626163     -0.069275]

if addprimary
	ross_secondary .+= ross_primary
end
