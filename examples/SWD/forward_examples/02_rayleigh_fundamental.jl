cd(@__DIR__)
ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"
using PythonCall
PD =  pyimport("disba").PhaseDispersion
np = pyimport("numpy")
thick = [10, 10, 15, 20, 20, 20, 20, 20, 0.]
vs = [3.38, 3.44, 3.66, 4.25, 4.35, 4.32, 4.315, 4.38, 4.5]
vpovervs = 1.77
vp = vs *vpovervs
ρ = 0.32 * vp .+ 0.77
periods = 10 .^range(log10(3), log10(80), 20)
## now obtain a forward dispersion 
pd = PD(np.array([thick vp vs ρ]).T...)
phase_vel = pd(np.array(periods), mode=0, wave="rayleigh").velocity
## plot
fig, ax = plt.subplots(1, 2, figsize=(8,4), gridspec_kw=Dict("width_ratios" => [1,2]))
ax[1].step(vs, [0, cumsum(thick)[1:end-1]...],)
ax[2].plot(periods, phase_vel, ".-")
ax[1].invert_yaxis()
ax[2].set_xlabel("Period s")
ax[2].set_ylabel("Phase velocity km/s")
ax[1].set_ylabel("Depth km")
ax[1].set_xlabel("Vs km/s")
ax[1].grid(); ax[2].grid()
fig.tight_layout()