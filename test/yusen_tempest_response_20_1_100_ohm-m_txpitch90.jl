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
tx_roll = 0.
tx_pitch = 90. #ACHTUNG!
tx_yaw = 0.
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
					  tx_yaw = tx_yaw)
# X and Z components in fT from GA-AEM
ross_primary = [50.4312      29.7847]
ross_secondary = [  -2.3613      -6.762
  -2.2662     -4.7924
  -2.1378     -4.0079
  -2.0113     -3.4924
  -1.8974     -3.1035
  -1.7945     -2.7646
  -1.6674     -2.3721
  -1.4912     -1.9079
  -1.2582     -1.4085
-0.97445    -0.92505
-0.67056    -0.52168
-0.40005    -0.24675
-0.20486   -0.097237
-0.091297   -0.032526
-0.034562   -0.008896]
