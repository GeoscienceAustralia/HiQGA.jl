using PyPlot, Revise, transD_GP, Random, Statistics, BenchmarkTools
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
zTx = -120
zRx = -80
x_rx = -115.
y_rx = 0.
rx_roll = 0.
rx_pitch = 0.
rx_yaw = 0.
tx_roll = 0.
tx_pitch = 90.
tx_yaw = 0.
##
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
## model and Ross' response - select one
#include("ross_tempest_response_10_ohm-m.jl")
#include("ross_tempest_response_100_ohm-m.jl")
#include("ross_tempest_response_20_10_1000_ohm-m.jl")
#include("yusen_tempest_response_20_1_100_ohm-m.jl")
include("yusen_tempest_response_20_1_100_ohm-m_txpitch90.jl");
## do it
transD_GP.TEMPEST1DInversion.getfieldTD!(tempest, zfixed, rho)
## plot
μ = transD_GP.AEM_VMD_HMD.μ
figure(figsize=(4,7))
s1 = subplot(211)
loglog(times, abs.(μ*tempest.Hz)*1e15, label="Bz", ".-")
loglog(times, abs.(μ*tempest.Hx)*1e15, label="Br", ".-")
xlabel("time s")
ylabel("B field 10⁻¹⁵ T")
ylim(1e-6, 50)
legend()
grid(true, which="both")
plot(times, abs.(ross_secondary[:,2]))
plot(times, abs.(ross_secondary[:,1]))
s2 = subplot(212, sharex=s1)
loglog(times, 100*abs.(abs.(μ*tempest.Hz*1e15)-abs.(ross_secondary[:,2]))./abs.(ross_secondary[:,2]),"x-")
loglog(times, 100*abs.(abs.(μ*tempest.Hx*1e15)-abs.(ross_secondary[:,1]))./abs.(ross_secondary[:,1]),"x-")
grid(true, which="both")
xlabel("time s")
ylabel("% diff")
plt.tight_layout()
## timing for many random layers
## time many random layers
doprofile=false
if doprofile
	using BenchmarkTools
	Random.seed!(435)
	nlayers = 50
	z = [-1e5, 0, 20, 50 .+ cumsum(15*rand(nlayers-3))...]
	ρ = [1e13, 10, 1, 10 .^(-0.5 .+ 1.5*rand(nlayers-3))...]
	@btime transD_GP.TEMPEST1DInversion.getfieldTD!($tempest, $z, $rho)
end
