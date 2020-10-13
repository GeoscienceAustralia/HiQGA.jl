srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, DelimitedFiles, Random, Statistics, Revise,
      AEM_VMD_HMD, SkyTEM1DInversion, GeophysOperator
## model fixed parts, i.e., air
Random.seed!(23)
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 100
# Note that the receiver and transmitter need to be in layer 1
zstart = 0.0
extendfrac, dz = 1.06, 2.
zall, znall, zboundaries = GeophysOperator.setupz(zstart, extendfrac, dz=dz, n=40)
z, ρ, nfixed = GeophysOperator.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
ρ[z.>=zstart] .= 10.
##  geometry and modeling parameters
using MAT
fdataname = "Line_115651_dbdt_gates_rangeidx_123.mat"
file = matopen(fdataname)
rRx = read(file, "rRx")
zRxLM = read(file, "LM_zRx")
zTxLM = read(file, "LM_zTx")
zRxHM = read(file, "HM_zRx")
zTxHM = read(file, "HM_zTx")
rTx = read(file, "rTxLoop")
lowpassfcs = read(file, "lowPassFilters")
ntimesperdecade = 10
nfreqsperdecade = 5
# Note that the receiver depth needs to be in same model layer as transmitter.
## LM times and ramp
LM_times = read(file, "LM_times")[:]
LM_ramp = read(file, "LM_ramp")
## HM times and ramp
HM_times = read(file, "HM_times")[:]
HM_ramp = read(file, "HM_ramp")
## LM operator
Flm = AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      lowpassfcs = lowpassfcs,
                      times  = LM_times,
                      ramp   = LM_ramp,
                      nmax   = nmax,
                      zTx    = zTxLM,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRxLM)
## HM operator
Fhm = AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      lowpassfcs = lowpassfcs,
                      times  = HM_times,
                      ramp   = HM_ramp,
                      nmax   = nmax,
                      zTx    = zTxHM,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRxHM)
## data and high altitude noise
LM_data = read(file, "d_LM")
HM_data = read(file, "d_HM")
LM_noise = read(file, "sd_LM")
HM_noise = read(file, "sd_HM")
## create operator
dlow, dhigh, σlow, σhigh = (LM_data, HM_data, LM_noise, HM_noise)./SkyTEM1DInversion.μ₀
aem = GeophysOperator.dBzdt(Flm, Fhm, vec(dlow), vec(dhigh),
                                  vec(σlow), vec(σhigh), z=z, ρ=ρ, nfixed=nfixed)
GeophysOperator.SkyTEM1DInversion.plotmodelfield!(Flm, Fhm, z, ρ, dlow, dhigh, σlow, σhigh;
                      figsize=(12,4), nfixed=nfixed, dz=dz, extendfrac=extendfrac)
