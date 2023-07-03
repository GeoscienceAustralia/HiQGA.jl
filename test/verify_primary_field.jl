using PyPlot, DelimitedFiles, Random, Statistics, HiQGA.transD_GP, Test

zfixed = [-1e5]
ρfixed = [1e12]

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
tx_yaw = -7
# electronics and stuff
include("tempest_electronics_halt.jl")
## create tempest operator 
tempest = transD_GP.TEMPEST1DInversion.Bfield(
    zTx = zTx, zRx = zRx, x_rx = x_rx, y_rx = y_rx,
    rx_roll = rx_roll, rx_pitch = rx_pitch, rx_yaw = rx_yaw,
    tx_roll = tx_roll, tx_pitch = tx_pitch, tx_yaw = tx_yaw,
	ramp = ramp, times = times,
	addprimary = true
)
geovec = transD_GP.TEMPEST1DInversion.extractnu(tempest)
transD_GP.TEMPEST1DInversion.returnprimary!(tempest, geovec)
μ = transD_GP.AEM_VMD_HMD.μ
@testset begin
    # values from Ross
    Bx, By, Bz = 30.2751, -3.7899, -13.8092
    @test isapprox(tempest.Hx[1]*μ*1e15, Bx, rtol=.001)
    @test isapprox(tempest.Hy[1]*μ*1e15, By, rtol=.001)
    @test isapprox(tempest.Hz[1]*μ*1e15, Bz, rtol=.001)
end     