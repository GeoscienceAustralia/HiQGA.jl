using PyPlot, DelimitedFiles, Random, Statistics, Test, HiQGA.transD_GP

Random.seed!(23)

zfixed = [-1e5]
ρfixed = [1e12]
nmax = 100

zstart = 0.
extendfrac, dz = 1.06, 2.
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=40, showplot=false)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
@info z

## geometry parameters for tempest
zTx = -120
zRx = -80
x_rx = -115.0
y_rx = 0.
rx_roll = 3.
rx_pitch = 6.
rx_yaw = 10.
tx_roll = 5.
tx_pitch = 4.
tx_yaw = -7.

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
ρ[1:end] .= ρ[1] # no contrasts
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
# plot before testing
transD_GP.TEMPEST1DInversion.plotmodelfield!(tempest,log10.(ρ[2:end]))

μ = transD_GP.TEMPEST1DInversion.μ₀
Bx, By, Bz = tempest.Hx[1]*μ*1e15, tempest.Hy[1]*μ*1e15, tempest.Hz[1]*μ*1e15
@testset "Primary field YPR order checks" begin
	@test isapprox(Bx, 30.2751, atol=1e-4)
	@test isapprox(By, -3.7899, atol=1e-4)
	@test isapprox(Bz, -13.8092, atol=1e-4)
end
