using Revise
cd(@__DIR__)
begin using Pkg; Pkg.activate("../../../") end
using PyPlot, DelimitedFiles, Random, Statistics, transD_GP

Random.seed!(23)

zfixed = [-1e5]
ρfixed = [1e12]
nmax = 100

zstart = 0.0
extendfrac, dz = 1.06, 2.
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=40)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
@info z

## geometry parameters for tempest
zTx = -120
zRx = -80
x_rx = -115.0
y_rx = 0.
rx_roll = 0.
rx_pitch = 0.
rx_yaw = 0.
tx_roll = 0.
tx_pitch = 0.
tx_yaw = 0.

#current ramp and gate times
ramp =  [  # -0.0400066666667    0.5
			-0.0400000000000    0.0
			-0.0399933333333   -0.5
			-0.0200066666667   -0.5
			-0.0200000000000    0.0
			-0.0199933333333    0.5
			-0.0000066666667    0.5
			 0.0000000000000    0.0
			 0.0000066666667   -0.5]

times = sqrt.([
0.0000066667*0.0000200000
0.0000333333*0.0000466667
0.0000600000*0.0000733333
0.0000866667*0.0001266667
0.0001400000*0.0002066667
0.0002200000*0.0003400000
0.0003533333*0.0005533333
0.0005666667*0.0008733333
0.0008866667*0.0013533333
0.0013666667*0.0021000000
0.0021133333*0.0032733333
0.0032866667*0.0051133333
0.0051266667*0.0079933333
0.0080066667*0.0123933333
0.0124066667*0.0199933333
])

## fill in detail in ohm-m
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50) .&(z.<80)] .= 1.
ρ[(z.>=80) .&(z.<100)] .= 20
ρ[(z.>=100) .&(z.<200)] .= 50
ρ[(z.>=200) .&(z.<250)] .= 80
ρ[(z.>=250)]            .= 150

## create total field operator (required for nuisance inversion)
tempest = transD_GP.TEMPEST1DInversion.Bfield(
    zTx = zTx, zRx = zRx, x_rx = x_rx, y_rx = y_rx,
    rx_roll = rx_roll, rx_pitch = rx_pitch, rx_yaw = rx_yaw,
    tx_roll = tx_roll, tx_pitch = tx_pitch, tx_yaw = tx_yaw,
	ramp = ramp, times = times,
	z=z,
	ρ=ρ,
	addprimary = true #this ensures that the geometry update actually changes everything that needs to be
)
# plot before adding noise
transD_GP.TEMPEST1DInversion.plotmodelfield!(tempest,z,ρ)
## compute noisy data to invert
transD_GP.TEMPEST1DInversion.set_noisy_data!(tempest, z, ρ)
