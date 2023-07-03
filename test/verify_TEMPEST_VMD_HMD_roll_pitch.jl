using PyPlot, Revise, HiQGA.transD_GP, Random, Statistics
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
## model and GA-AEM TEMPEST responses
# GA-AEM does not model yaw as far as I can see so leaving those examples out
# TX and RX pitch are the big hitters anyway
fnames=[
"yusen_tempest_response_20_1_100_ohm-m_pure_VMD.jl",
"yusen_tempest_response_20_1_100_ohm-m_txpitch90.jl",
"yusen_tempest_response_20_1_100_ohm-m_txroll10_txpitch10.jl",
"yusen_tempest_response_20_1_100_ohm-m_rxroll10_rxpitch10.jl",
#"yusen_tempest_response_20_1_100_ohm-m_rxyaw10.jl"
#"yusen_tempest_response_20_1_100_ohm-m_rxroll10_rxpitch10_rxyaw10.jl"
]
## model conductivities and yaw pitch roll in file
for (ifn, fn) in enumerate(fnames)
	include(fn)
	transD_GP.TEMPEST1DInversion.getfieldTD!(tempest, zfixed, rho)
	μ = transD_GP.AEM_VMD_HMD.μ
	figure(figsize=(8,9))
	s = subplot(121)
	s.step(vcat(rho[2:end],rho[end]), vcat(zfixed[2:end], zfixed[end]+50))
	s.set_xscale("log")
	s.invert_xaxis()
	xlabel("ρ")
	ylabel("depth m")
	grid(true, which="both")
	s1 = subplot(222)
	loglog(times, abs.(μ*tempest.Hz)*1e15, label="Bz", ".-")
	loglog(times, abs.(μ*tempest.Hx)*1e15, label="Bx", ".-")
	xlabel("time s")
	ylabel("B field 10⁻¹⁵ T")
	# ylim(1e-6, 50)
	legend()
	grid(true, which="both")
	plot(times, abs.(ross_secondary[:,2]))
	plot(times, abs.(ross_secondary[:,1]))
	s2 = subplot(224, sharex=s1)
	loglog(times, 100*abs.(μ*tempest.Hz*1e15-ross_secondary[:,2])./abs.(ross_secondary[:,2]),"x-")
	loglog(times, 100*abs.(μ*tempest.Hx*1e15-ross_secondary[:,1])./abs.(ross_secondary[:,1]),"x-")
	grid(true, which="both")
	xlabel("time s")
	ylabel("% diff")
	plt.suptitle("$fn:\ntx_pitch:$tx_pitch, tx_roll:$tx_roll,  tx_yaw=$tx_yaw,\nrx_pitch:$rx_pitch, rx_roll:$rx_roll,  rx_yaw=$rx_yaw")
	plt.tight_layout()
end
