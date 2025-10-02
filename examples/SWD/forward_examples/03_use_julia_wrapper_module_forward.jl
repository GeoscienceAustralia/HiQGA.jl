cd(@__DIR__) # Change to this directory
ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"
include("../src/SWDPhysics.jl")
np = SWDPhysics.np
pyconvert = SWDPhysics.pyconvert

## model
thick = [10., 10., 10, 10,  10,  10, 10,  10, 10]
vp        = [7.0 , 6.8, 7., 7.6, 8.4, 9., 9.4, 9.6, 9.5]  
vs        = [3.5,  3.4, 3.5,  3.8,  4.2,  4.5,  4.7,  4.8,  4.75]
ρ         = [2., 2., 2., 2., 2., 2., 2., 2., 2.]

## convert to numpy arrays
periodswanted = np.array(10 .^range(0.0, 3.0, 100))
modeswanted = np.array([0, 1, 2])
waveswanted = pyconvert.(String, (["rayleigh", "love"]))
# one array for each computation call
periods = reduce(vcat, [[periodswanted for modes in modeswanted] for wave in waveswanted])
modes = reduce(vcat, [[mode for mode in modeswanted] for wave in waveswanted])
waves = reduce(vcat, [[wave for mode in modeswanted] for wave in waveswanted])

## compute
swdresults = SWDPhysics.getdispersion(thick, vp, vs, ρ, periods, modes, waves);

## plot
using PyPlot
fig, ax = plt.subplots(1, 2, sharex=true, sharey=true, figsize=(8,4))
for (s) in swdresults
    i = s.wave == "rayleigh" ? 1 : 2
    ax[i].plot(s.period, s.velocity, label="mode=$(s.mode+1)")
end
ax[1].set_xscale("log")
ax[1].set_ylabel("Vₚ km/s")
[(ax[i].set_xlabel("Period s"); ax[i].grid(axis="y"); ax[i].legend(fancybox=true, framealpha=0.5)) for i in 1:2]
ax[1].set_title("Rayleigh")
ax[2].set_title("Love")
fig.tight_layout()
