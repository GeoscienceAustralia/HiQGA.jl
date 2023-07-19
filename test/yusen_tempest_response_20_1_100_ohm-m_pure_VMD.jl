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
rx_roll = 0. #ACHTUNG!
rx_pitch = 0. #ACHTUNG!
rx_yaw = 0. #ACHTUNG!
tx_roll = 0. #ACHTUNG!
tx_pitch = 0. #ACHTUNG!
tx_yaw = 0. #ACHTUNG!
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
					  addprimary = addprimary,
					  z = zfixed,
					  ρ = rho)
# Fields in fT from GA-AEM
ross_primary = [25.7847     -18.7313]
ross_secondary = [   6.762     -8.0374
  4.7924     -6.6042
  4.0079     -5.9113
  3.4924     -5.3952
  3.1035     -4.9777
  2.7646      -4.608
  2.3721     -4.1681
  1.9079     -3.6087
  1.4085     -2.9386
0.92505     -2.1946
0.52168     -1.4595
0.24675    -0.84605
0.097237    -0.42402
0.032526    -0.18624
0.008896   -0.069836]

if addprimary
	ross_secondary .+= ross_primary
end
