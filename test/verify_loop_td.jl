using PyPlot, Revise, transD_GP, Random, SpecialFunctions
## model
zfixed   = [-1e5,   0,    ]
rho      = [1e20,   1,  ]
nmax = 200
##  geometry
rTx = 15
rRx = 0.001
zRx = -0.001
zTx = -0.002
doconvramp = false
modelprimary = true
nkᵣeval = 200
ntimesperdecade = 15
nfreqsperdecade = 15
times = 10 .^LinRange(-6,1,40)
## build operator
F = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      times  = times,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx,
                      doconvramp = doconvramp,
                      modelprimary = modelprimary,
                      nkᵣeval = nkᵣeval)
## get TD field
transD_GP.AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
## Analytic loop
μ = transD_GP.AEM_VMD_HMD.μ
sig0 = 1/rho[2]
theta = sqrt.(μ*sig0./4F.interptimes)
a     = rTx
dbzdt = -1 ./(sig0*a^3)*(3*erf.(theta*a) - 2/sqrt(pi)*theta.*a.*(3 .+ 2*theta.^2*a^2).*exp.(-theta.^2*a^2))
dbzdt = dbzdt/(pi*a^2)
##
figure()
s1 = subplot(211)
loglog(F.interptimes,μ*abs.(F.HTD_z_interp), label="loop W&H")
loglog(F.interptimes,abs.(dbzdt), "--", label="loop W&H transient half space ")
xlim(extrema(times))
ylabel("dBzdt V/Am^4")
grid(which="both")
legend()
s2 = subplot(212, sharex=s1)
# loglog(F.interptimes, 100*abs.(-AEM_VMD_HMD.μ*F.HTDinterp-dbzdt)./abs.(dbzdt))
# signs are supposed to be the same at loop centre ...
loglog(F.interptimes, 100*abs.(abs.(μ*F.HTD_z_interp)-abs.(dbzdt))./abs.(dbzdt))
grid(which="both")
ylabel("% difference")
xlabel("time s")
