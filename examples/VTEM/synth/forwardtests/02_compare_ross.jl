using PyPlot, DelimitedFiles, Random, Statistics,
      HiQGA.transD_GP
## model fixed parts, i.e., air
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 200
# Note that the receiver and transmitter need to be in layer 1
zstart = 0.0
extendfrac, dz = 1.03, 1.5
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=65, showplot=false)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
##  geometry and modeling parameters
rRx = 0.0
zRx = -30.0
zTx = -30.0
freqlow = 1e-4
freqhigh = 1e6
nkᵣeval = 50
ntimesperdecade = 10
nfreqsperdecade = 5
calcjacobian = true
include("../waveletapprox/electronics_halt.jl")
## LM operator
F = transD_GP.AEM_VMD_HMD.HFieldDHT(;
                      times,
                      ramp,
                      nmax,
                      zTx,
                      rRx,
                      rTx,
                      zRx,
                      freqlow,
                      freqhigh,
                      nkᵣeval,
                      ntimesperdecade,
                      nfreqsperdecade,
                      calcjacobian
                      )
# Ross responses
fnames = ["ross_vtem_response_100_ohm-m.jl",
          "ross_vtem_response_100_1_1000_ohm-m.jl"]
# get field
μ = transD_GP.AEM_VMD_HMD.μ
##
for fn in fnames
    include(fn)
    @time transD_GP.AEM_VMD_HMD.getfieldTD!(F, zfixed, ρ)
    figure(figsize=(6,7))
	s = subplot(121)
	s.step(vcat(ρ[2:end],ρ[end]), vcat(zfixed[2:end], zfixed[end]+50))
	s.set_xscale("log")
	s.invert_xaxis()
    s.invert_yaxis()
	xlabel("ρ")
	ylabel("depth m")
	grid(true, which="both")
	s1 = subplot(222)
	loglog(F.times, abs.(μ*F.dBzdt)*1e12, label="Bz Julia", ".-")
	xlabel("time s")
	ylabel("dB/dt field pV/Am^4")
	# ylim(1e-6, 50)
	plot(times, abs.(ross_secondary[:,2]), label="GA-AEM")
	legend()
	grid(true, which="both")
	s2 = subplot(224, sharex=s1)
	loglog(F.times, 100*abs.(μ*F.dBzdt*1e12+ross_secondary[:,2])./abs.(ross_secondary[:,2]),"x-")
	grid(true, which="both")
	xlabel("time s")
	ylabel("% diff")
	plt.suptitle("Pure VMD no roll or pitch")
	plt.draw()
	plt.tight_layout()
end    