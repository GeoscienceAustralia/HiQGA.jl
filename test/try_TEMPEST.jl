using PyPlot, Revise, transD_GP, Random, Statistics
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
rRx = 115.0
modelprimary = false
getradialH = true
provideddt = false
Ftempest = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
					  nkᵣeval = nkᵣeval,
                      times  = times,
                      ramp   = ramp,
                      zTx    = zTx,
                      rRx    = rRx,
                      zRx    = zRx,
					  modelprimary = modelprimary,
					  getradialH = getradialH,
					  provideddt = provideddt)
## model and Ross' response - select one
include("ross_tempest_response_10_ohm-m.jl")
#include("ross_tempest_response_100_ohm-m.jl")
#include("ross_tempest_response_20_10_1000_ohm-m.jl")
## do it
transD_GP.AEM_VMD_HMD.getfieldTD!(Ftempest, zfixed, rho)
## plot
μ = transD_GP.AEM_VMD_HMD.μ
figure(figsize=(4,7))
s1 = subplot(211)
loglog(Ftempest.times, abs.(μ*Ftempest.dBzdt)*1e15, label="Bz", ".-")
loglog(Ftempest.times, abs.(μ*Ftempest.dBrdt)*1e15, label="Br", ".-")
xlabel("time s")
ylabel("B field 10⁻¹⁵ T")
ylim(1e-6, 50)
legend()
grid(true, which="both")
plot(times, abs.(ross_secondary[:,2]))
plot(times, abs.(ross_secondary[:,1]))
s2 = subplot(212, sharex=s1)
loglog(times, 100*abs.(μ*Ftempest.dBzdt*1e15-ross_secondary[:,2])./abs.(ross_secondary[:,2]),"x-")
loglog(times, 100*abs.(μ*Ftempest.dBrdt*1e15-ross_secondary[:,1])./abs.(ross_secondary[:,1]),"x-")
grid(true, which="both")
xlabel("time s")
ylabel("% diff")
plt.tight_layout()
