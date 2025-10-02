module SWDPhysics
using PythonCall
PD =  pyimport("disba").PhaseDispersion
np = pyimport("numpy")

mutable struct SWDispersion
    period  :: Array{Float64, 1}
    velocity :: Array{Float64, 1}
    mode     :: Int
    wave     :: String
    type     :: String
end

## convert to julia
function cp_ptoj(cp)
    period   = pyconvert(Array{Float64, 1}, cp.period)
    velocity = pyconvert(Array{Float64, 1}, cp.velocity)
    mode     = pyconvert(Int, cp.mode)
    wave     = pyconvert(String, cp.wave)
    type     = pyconvert(String, cp.type)
    SWDispersion(period, velocity, mode, wave, type)
end

function getdispersion(thick, vp, vs, ρ, 
    periods, # array of periods, as many arrays as dispersion curves
    modes, # as many as dispersion curves
    waves # as many as dispersion curves
)
    # now obtain a forward dispersion object
    pd = PD(np.array([thick vp vs ρ]).T...)
    map(zip(periods, modes, waves)) do (t, mode, wave)
        cp_ptoj(pd(t, mode=mode, wave=wave))
    end
end

end