cd(@__DIR__) # Change to this directory
ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"
include("../src/SWDPhysics.jl")
include("../src/SWDInversion.jl")
using .SWDInversion, HiQGA, NearestNeighbors, PyPlot, Random
Random.seed!(11)
np = SWDPhysics.np
pyconvert = SWDPhysics.pyconvert

## model 
zgiven  = cumsum([10, 10, 15, 20, 20, 20, 20, 20, 0.])
vsgiven = [3.38, 3.44, 3.66, 4.25, 4.35, 4.32, 4.315, 4.38, 4.5]
# convert to increasing fine grid depth form
zstart, extendfrac, dz = 0., 1.05, 0.75
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50)
thick = [diff(zboundaries)..., 0.]
vs = transD_GP.gridpoints(zgiven, zall, vsgiven)
fig, ax = plt.subplots(figsize=(3,4))
ax.step(vs, zboundaries)
ax.invert_yaxis(); ax.set_ylabel("Depth km"); ax.set_xlabel("Vs km/s")
fig.tight_layout()
#
vpovervs = 1.77
vp = vs * vpovervs
ρ = 0.32 * vp .+ 0.77
σ =  0.02 # noise in km/s
## convert to numpy arrays
periodswanted = np.array(10 .^range(log10(3), log10(80), 20))
modeswanted = np.array([0])
waveswanted = pyconvert.(String, (["rayleigh"]))
# one array for each computation call
periods = reduce(vcat, [[periodswanted for modes in modeswanted] for wave in waveswanted])
modes = reduce(vcat, [[mode for mode in modeswanted] for wave in waveswanted])
waves = reduce(vcat, [[wave for mode in modeswanted] for wave in waveswanted])

## compute and add noise
swdresults = SWDPhysics.getdispersion(thick, vp, vs, ρ, periods, modes, waves);
noisy_data = [s.velocity + σ*randn(size(s.velocity)) for s in swdresults]
noisesd = [σ*ones(size(s.velocity)) for s in swdresults]

## plot
using PyPlot
fig, ax = plt.subplots(1, 2, sharex=true, sharey=true, figsize=(8,4))
for (s, nd) in zip(swdresults, noisy_data)
    i = s.wave == "rayleigh" ? 1 : 2
    ax[i].plot(s.period, s.velocity, label="mode=$(s.mode+1)")
    ax[i].plot(s.period, nd, ".")
end
ax[1].set_ylabel("Vₚ km/s")
[(ax[i].set_xlabel("Period s"); ax[i].grid(axis="y"); ax[i].legend(fancybox=true, framealpha=0.5)) for i in 1:2]
ax[1].set_title("Rayleigh")
ax[2].set_title("Love")
fig.tight_layout()

## make a SWD struct
swd = SWDInversion.SWDispInversion(;
    thick, vpovervs, ρ,     
    d = noisy_data, σ = noisesd,
    periods, modes, waves
)
## check chi^2 error
r = reduce(vcat, swd.d) - reduce(vcat, getfield.(swdresults, :velocity))
n = reduce(vcat, swd.σ)
@info "χ² is $(sum((r./n).^2)) for $(length(r)) data"