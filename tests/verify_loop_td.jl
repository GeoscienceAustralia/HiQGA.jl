srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Revise, AEM_VMD_HMD, Random, SpecialFunctions
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
ntimesperdecade = 15
nfreqsperdecade = 15
times = 10 .^LinRange(-6,1,40)
## build operator
F = AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      times  = times,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx,
                      doconvramp = doconvramp)
## get TD field
AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
## Analytic loop
sig0 = 1/rho[2]
theta = sqrt.(4*pi*1e-7*sig0./4F.interptimes)
a     = rTx
dbzdt = -1 ./(sig0*a^3)*(3*erf.(theta*a) - 2/sqrt(pi)*theta.*a.*(3 .+ 2*theta.^2*a^2).*exp.(-theta.^2*a^2))
dbzdt = dbzdt/(pi*a^2)
##
figure()
loglog(F.interptimes,4*pi*1e-7*abs.(F.HTDinterp), label="circle W&H")
loglog(F.interptimes,abs.(dbzdt), "--", label="circle Christiansen")
xlim(extrema(times))
grid()
legend()
