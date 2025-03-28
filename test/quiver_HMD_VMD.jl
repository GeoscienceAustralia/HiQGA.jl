using PyPlot, HiQGA.transD_GP, Random, Statistics
## waveform and times
ramp =  [  # -0.0400066666667    0.5
			-0.0400000000000    0.0
			-0.0399933333333   -0.5
			-0.0200066666667   -0.5
			-0.0200000000000    0.0
			-0.0199933333333    0.5
			-0.0000066666667    0.5
			 0.0000000000000    0.0
			 0.0000066666667   -0.5]

times = vec(10 .^mean(log10.([
            0.0000066667	0.0000200000
			0.0000333333	0.0000466667
			0.0000600000	0.0000733333
			0.0000866667	0.0001266667
			0.0001400000	0.0002066667
			0.0002200000	0.0003400000
			0.0003533333	0.0005533333
			0.0005666667	0.0008733333
			0.0008866667	0.0013533333
			0.0013666667	0.0021000000
			0.0021133333	0.0032733333
			0.0032866667	0.0051133333
			0.0051266667	0.0079933333
			0.0080066667	0.0123933333
			0.0124066667	0.0199933333]), dims=2))
## Set up operator
ntimesperdecade = 10
nfreqsperdecade = 5
nkᵣeval = 60
# x,y,z in z downward!!
# +x is flight direction
zTx = -120
zRx = -80
x_rx = -115.
y_rx = 0.
# these are in z upward, ga-aem convention
rx_roll = 0.
rx_pitch = 0.
rx_yaw = 0.
tx_roll = 0. # GA AEM downward pointing dipole is roll 180
tx_pitch = 90.
tx_yaw = 0.
# model
zfixed   = [-1e6,   0, 40, 60]
rho      = [1e12,   20, 1, 100]
## do it
xmin, xmax, dx = -4.25, 4.25, 0.25
ymin, ymax, dy = -3.25, 3.25, 0.25
xrx = xmin:dx:xmax
yrx = ymin:dy:ymax
u = zeros(length(yrx), length(xrx), length(times))
v = zeros(length(yrx), length(xrx), length(times))
w = zeros(length(yrx), length(xrx), length(times))
for (jx, x) in enumerate(xrx), (iy, y) in enumerate(yrx)
	tempest = transD_GP.TEMPEST1DInversion.Bfield(
  					  ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
					  nkᵣeval = nkᵣeval,
                      times  = times,
                      ramp   = ramp,
					  zTx    = zTx,
					  zRx    = zRx,
					  x_rx   = x,
					  y_rx   = y,
					  rx_roll = rx_roll,
					  rx_pitch = rx_pitch,
					  rx_yaw = rx_yaw,
					  tx_roll = tx_roll,
					  tx_pitch = tx_pitch,
					  tx_yaw = tx_yaw,
					  strictgeometry = false)
	transD_GP.TEMPEST1DInversion.getfieldTD!(tempest, zfixed, rho)
	# Roll180 underneath to get back to Z down co-ordinate system, easier to plot
	# and make sense
	u[iy,jx,:] = tempest.Hx
	v[iy,jx,:] = -tempest.Hy
	w[iy,jx,:] = -tempest.Hz
end
## plot it
f, ax = plt.subplots(1,3, sharex=true, sharey=true, figsize=(16,4))

sleeptime = 1
ntimes = 15
for it = 1:ntimes
	im1 = ax[1].pcolormesh(xrx,yrx, u[:,:,it])
	cb1 = colorbar(im1,ax=ax[1])
	ax[1].quiver(xrx, yrx, u[:,:,it], zeros(size(u,1), size(u,2)), width=0.005)
	ax[1].set_xlabel("x m")
	ax[1].set_ylabel("y m")
	ax[1].set_title("Hx")
	im2 = ax[2].pcolormesh(xrx,yrx, v[:,:,it])
	cb2 = colorbar(im2,ax=ax[2])
	ax[2].quiver(xrx, yrx, zeros(size(v,1), size(v,2)), v[:,:,it], width=0.005)
	ax[2].set_xlabel("x m")
	ax[2].set_ylabel("y m")
	ax[2].set_title("Hy")
	im3 = ax[3].pcolormesh(xrx,yrx, w[:,:,it])
	cb3 = colorbar(im3,ax=ax[3])
	ax[3].set_xlabel("x m")
	ax[3].set_ylabel("y m")
	ax[3].set_title("Hz")
	plt.draw()
	sleep(sleeptime)
	if it != ntimes
		ax[1].clear()
		ax[2].clear()
		ax[3].clear()
		cb1.remove()
		cb2.remove()
		cb3.remove()
	end
end
