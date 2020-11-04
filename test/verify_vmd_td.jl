srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Revise, AEM_VMD_HMD, Random
## set up
nfreqsperdecade     = 20
ntimesperdecade     = 20
modelprimary        = true
provideddt          = true
doconvramp          = false
nkᵣeval             = 200
times               = 10 .^LinRange(-5,-1.1, 60)
lowpassfcs          = [1e6]
freqhigh = 1e6
freqlow = 1e-8
## model
zfixed   = [-1e5,   0]
rho      = [1e12,   100]
nmax = 200
##  geometry
rRx = 100.
zRx = 0.
zTx = 0.
modelprimary = true
##
F = AEM_VMD_HMD.HFieldDHT(freqlow = freqlow,
                      freqhigh = freqhigh,
                      zTx    = zTx,
                      rRx    = rRx,
                      times = times,
                      nfreqsperdecade = nfreqsperdecade,
                      ntimesperdecade = ntimesperdecade,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary,
                      provideddt = provideddt,
                      lowpassfcs = lowpassfcs)
AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
## plot FD
figure()
loglog(F.freqs, abs.(real.(F.HFD_z)))
loglog(F.freqs, abs.(imag.(F.HFD_z)))
xlim(1e-1,1e5)
ylim(1e-12, 1e-6)
xlabel("Hz")
title("Compare with W&H 4.2")
## plot time domain
figure(figsize=(4,7))
# first dhzdt
F.provideddt = true
AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
loglog(F.interptimes*1e3, abs.(F.HTD_z_interp), label="dhz/dt")
xlim(1e-5,1e3)
ylim(1e-11,1e-1)
grid(true, which="both")
# then dhdt
F.provideddt = false
AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
loglog(F.interptimes*1e3, abs.(F.HTD_z_interp), label="hz")
legend()
xlabel("time in ms")
title("Compare with W&H Fig 4.4")
plt.tight_layout()
## Radial fields
figure(figsize=(4,7))
F.useprimary = 0.
F.getradialH = true
F.provideddt = true
AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
loglog(F.interptimes*1e3, abs.(F.HTD_r_interp), label="dhr/dt")
ylim(1e-11,1e-1)
legend()
ylabel("dh/dt")
# grid(true, which="both")
xlabel("time in ms")
F.provideddt = false
AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
ax = twinx(gca())
ax.loglog(F.interptimes*1e3, abs.(F.HTD_r_interp), label="hr", color="k")
ax.set_ylabel("h")
ax.set_xlim(1e-3,10^3.5)
ax.set_ylim(10^-7.05,10^-16.9)
ax.invert_yaxis()
# ax.grid(true, which="both")
title("Compare with W&H Fig 4.5")
legend(loc="lower left")
plt.tight_layout()
## investigative function
function plotkernelfunc(G::typeof(F); decfactor=10, J0=true)
    figure()
    for ifreq = 1:decfactor:length(G.freqs)
        freq = G.freqs[ifreq]
        if J0
            loglog(G.interpkᵣ, abs.(real(G.J0_kernel_v[:,ifreq])), label="$freq Hz", "--")
            loglog(G.interpkᵣ, abs.(imag(G.J0_kernel_v[:,ifreq])), "-")
        else
            loglog(G.interpkᵣ, abs.(real(G.J1_kernel_v[:,ifreq])), label="$freq Hz", "--")
            loglog(G.interpkᵣ, abs.(imag(G.J1_kernel_v[:,ifreq])), "-")
        end
    end
    J0 ? title("J0 v") : title("J1 v")
    legend()
    figure()
    s1 = subplot(211)
    if J0
        pcolormesh(G.interpkᵣ, G.freqs, log10.(abs.(real(G.J0_kernel_v)))')
    else
        pcolormesh(G.interpkᵣ, G.freqs, log10.(abs.(real(G.J1_kernel_v)))')
    end
    colorbar()
    J0 ? title("real J0 kernel v") : title("real J1 kernel v")
    ylabel("Hz")
    subplot(212, sharex=s1, sharey=s1)
    if J0
        pcolormesh(G.interpkᵣ, G.freqs, log10.(abs.(imag(G.J0_kernel_v)))')
    else
        pcolormesh(G.interpkᵣ, G.freqs, log10.(abs.(imag(G.J1_kernel_v)))')
    end
    colorbar()
    J0 ? title("imag J0 kernel v") : title("imag J1 kernel v")
    gca().set_yscale("log")
    gca().set_xscale("log")
    xlabel("kᵣ 1/m"); ylabel("Hz")
    plt.tight_layout()
end
