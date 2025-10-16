cd(@__DIR__)
ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"
using PythonCall
PD =  pyimport("disba").PhaseDispersion
np = pyimport("numpy")

# Velocity model
# thickness, Vp, Vs, density
# km, km/s, km/s, g/cm3
velocity_model = np.array([
   [10.0, 7.00, 3.50, 2.00],
   [10.0, 6.80, 3.40, 2.00],
   [10.0, 7.00, 3.50, 2.00],
   [10.0, 7.60, 3.80, 2.00],
   [10.0, 8.40, 4.20, 2.00],
   [10.0, 9.00, 4.50, 2.00],
   [10.0, 9.40, 4.70, 2.00],
   [10.0, 9.60, 4.80, 2.00],
   [10.0, 9.50, 4.75, 2.00],
])

# Periods must be sorted starting with low periods
t = np.logspace(0.0, 3.0, 100)

# Compute the 3 first Rayleigh- and Love- wave modal dispersion curves
# Fundamental mode corresponds to mode 0
pd = PD(np.array(pylist(velocity_model)).T...) # change *velocity_model.T using list and splat...
cpr = [pd(t, mode=i, wave="rayleigh") for i in pyrange(3)] # range to pyrange
cpl = [pd(t, mode=i, wave="love") for i in pyrange(3)] # range to pyrange

## convert to julia
function cp_ptoj(cp_py)
    map(cp_py) do cp
        (period   = pyconvert(Array{Float64}, cp.period), 
         velocity = pyconvert(Array{Float64}, cp.velocity),
         mode     = pyconvert(Int, cp.mode),
         wave     = pyconvert(String, cp.wave),
         type     = pyconvert(String, cp.type))
    end
end

cpr_j, cpl_j = cp_ptoj.([cpr, cpl])

## plot
using PyPlot
fig, ax = plt.subplots(1, 2, sharex=true, sharey=true, figsize=(8,4))
for (r,l) in zip(cpr_j, cpl_j)
    ax[1].plot(r.period, r.velocity, label="mode=$(r.mode+1)")
    ax[2].plot(l.period, l.velocity, label="mode=$(l.mode+1)")
end
ax[1].set_xscale("log")
ax[1].set_ylabel("Vₚ km/s")
[(ax[i].set_xlabel("Period s"); ax[i].grid(axis="y"); ax[i].legend(fancybox=true, framealpha=0.5)) for i in 1:2]
ax[1].set_title("Rayleigh")
ax[2].set_title("Love")
fig.tight_layout()
## 