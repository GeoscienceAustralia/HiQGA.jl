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
rx_roll = 10. #ACHTUNG!
rx_pitch = 10. #ACHTUNG!
rx_yaw = 10. #ACHTUNG!
tx_roll = 0.
tx_pitch = 0.
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
					  addprimary=addprimary)
# X and Z components in fT from GA-AEM
ross_primary = [28.5962      -13.689]
ross_secondary = [
   8.03371     -6.620856
  5.849011     -5.572887
  4.957888     -5.037114
  4.362005     -4.626088
  3.907595     -4.288714
  3.510647     -3.988984
  3.048846      -3.63051
  2.496073     -3.168608
  1.889611     -2.605418
  1.286286     -1.967748
0.7633423     -1.324871
0.3876874    -0.7776925
0.1682725    -0.3943537
0.06388109    -0.1749779
0.02070349   -0.06618546]

if addprimary
	ross_secondary .+= ross_primary
end
