using PyPlot, Revise, transD_GP, Random, Statistics
include("skytem_response.jl")
##  modeling parameters
ntimesperdecade = 10
nfreqsperdecade = 5
## LM operator
Flm = transD_GP.SkyTEM1DInversion.AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      lowpassfcs = lowpassfcs,
                      times  = LM_times,
                      ramp   = LM_ramp,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRxLM)
## HM operator
Fhm = transD_GP.SkyTEM1DInversion.AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      lowpassfcs = lowpassfcs,
                      times  = HM_times,
                      ramp   = HM_ramp,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRxHM)
## get the fields in place
@time transD_GP.SkyTEM1DInversion.AEM_VMD_HMD.getfieldTD!(Flm, zfixed, rho)
@time transD_GP.SkyTEM1DInversion.AEM_VMD_HMD.getfieldTD!(Fhm, zfixed, rho)
## plot
figure()
s1 = subplot(211)
loglog(Flm.times,4*pi*1e-7*abs.(Flm.dBzdt), label="lm_julia")
loglog(Fhm.times,4*pi*1e-7*abs.(Fhm.dBzdt), label="hm_julia")
loglog(Flm.times,dLM, ".", label="lm_skytem")
loglog(Fhm.times,dHM, ".", label="hm_skytem")
legend()
grid(which="both")
ylabel("dBzdt V/Am^4")
s2 = subplot(212, sharex=s1)
semilogx(Flm.times,(4*pi*1e-7*Flm.dBzdt-dLM)./dLM*100, label="lm")
semilogx(Fhm.times,(4*pi*1e-7*Fhm.dBzdt-dHM)./dHM*100, label="hm")
ylabel("% difference")
xlabel("time s")
legend()
grid(which="both")
