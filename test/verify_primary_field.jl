using PyPlot, DelimitedFiles, Random, Statistics, HiQGA.transD_GP, Test

zfixed = [-1e5]
ρfixed = [1e12]

## geometry parameters for tempest 
# z down position coordinates
# GA-AEM z up rotation coordinates
# +x is flight direction
# this is a bad geovec, not true, updated version underneath
zTx      = -110
zRx      = -70
x_rx     = -113.0
y_rx     = 2.
rx_roll  = 4.
rx_pitch = 7.
rx_yaw   = 8.
tx_roll  = 7.
tx_pitch = 3.
tx_yaw   = -2
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
# true geovec from Ross
geovec = [-120, -80, -115, 0, 3, 6, 10, 5, 4, -7.]
transD_GP.TEMPEST1DInversion.returnprimary!(tempest, geovec)
μ = transD_GP.AEM_VMD_HMD.μ
@testset begin
    # values from Ross in fT
    Bx, By, Bz = 30.2751, -3.7899, -13.8092
    @test isapprox(tempest.Hx[1]*μ*1e15, Bx, rtol=.001)
    @test isapprox(tempest.Hy[1]*μ*1e15, By, rtol=.001)
    @test isapprox(tempest.Hz[1]*μ*1e15, Bz, rtol=.001)
end     