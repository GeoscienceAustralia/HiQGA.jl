srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Revise, AEM_VMD_HMD, Random, Statistics
## waveform and times
ramp = [    -0.0200000000000    0.0
			-0.0199933333333    1.0
			-0.0000066666667    1.0
			 0.0000000000000    0.0]
			 # 0.0000066666667   -1.0
			 # 0.0199933333333   -1.0
			 # 0.0200000000000    0.0]

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
nkᵣeval = 50
zTx = -120
zRx = -80
rRx = 115.0
modelprimary = false
getradialH = true
provideddt = false
Ftempest = AEM_VMD_HMD.HFieldDHT(
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
## model
zfixed   = [-1e6,   0, 400]
rho      = [1e12,   100, 100]
## do it
AEM_VMD_HMD.getfieldTD!(Ftempest, zfixed, rho)
## plot
figure(figsize=(4,7))
loglog(Ftempest.times, abs.(AEM_VMD_HMD.μ*Ftempest.dBzdt)*1e15, label="Bz", ".-")
loglog(Ftempest.times, abs.(AEM_VMD_HMD.μ*Ftempest.dBrdt)*1e15, label="Br", ".-")
ylim(1e-4, 10)
legend()
grid(true, which="both")
