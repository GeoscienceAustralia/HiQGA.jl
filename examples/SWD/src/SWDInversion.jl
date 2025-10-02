# make sure to include this file 
# AFTER SWDPhysics.jl has been included
module SWDInversion
using HiQGA.transD_GP.AbstractOperator
# using AbstractOperator (above) says that transD_GP can 
# use this kind of physics, through a misfit function

import HiQGA.transD_GP.AbstractOperator.get_misfit
# get_misfit (above) extends all the other misfit functions 
# in transD_GP. Because we are extending a method
# we must use import instead of using.

import HiQGA.transD_GP.Model, HiQGA.transD_GP.Options
# the above line now recognizes GP models and McMC options defined in src/TransD_GP_MCMC.jl

# The above three lines are the bare minimum we need to use transD_GP

# But wait! we need to specialize this inversion module 
# for BarPhysics ... so let's do that. This is basically the 
# geophysics "forward call"
using ..SWDPhysics  
import ..SWDPhysics.getdispersion
np = SWDPhysics.np

# Now let's actually define our SWDInversion operator and for simplicity
# say it is a subtype of 1D operator (operators are defined in src/AbstractOperator.jl)
mutable struct SWDispInversion<:Operator1D
    # thickness is one fixed vector
    thick :: Array{Float64, 1}
    vpovervs :: Float64
    ρ :: Array{Float64, 1}
    # all these are vector of vector
    d :: Vector{Array{Float64, 1 }} # observations at different periods for each mode
    σ :: Vector{Array{Float64, 1 }} # noise for different periods at each mode
    r
    periods 
    modes 
    waves
end   

function SWDispInversion(;thick=nothing, vpovervs=nothing, ρ=nothing, d=nothing, σ=nothing, 
                        periods=nothing, modes=nothing, waves=nothing)
    for x in (thick, vpovervs, ρ, d, σ, periods, modes, waves)
        @assert !isnothing(x)
    end
    r = [zeros(length(p)) for p in periods]
    SWDispInversion(thick, vpovervs, ρ, d, σ, r, periods, modes, waves)
end
# Now let's define (extend) a misfit function for
# a 1D model parameter vector m, the physics for which 
# can be computed through methods in SWDPhysics
# using an operator of type SWDispersion, using McMC prior bounds 
# and options stored in a struct of type Options
function get_misfit(m::Model, opt::Options, swd::SWDispInversion)
    get_misfit(m.fstar, opt, swd) # dispatches to helper function below
end
# above defined function and type signature MUST be defined

function get_misfit(x::AbstractArray, opt::Options, swd::SWDispInversion)
    # m.fstar is the 1D array interpolated by the GP mean
    # at locations opt.xall
    chi2by2 = 0.
    opt.debug && (return chi2by2)
    chi2by2 = get_misfit(x, swd)
end

function get_misfit(x::AbstractArray, swd::SWDispInversion)
    vs = vec(x) # nasty
    try
        disps = SWDPhysics.getdispersion(swd.thick, vs*swd.vpovervs, vs, swd.ρ, swd.periods, swd.modes, swd.waves)
        map(zip(disps, swd.d, swd.σ, swd.r)) do (disp, d, sd, res)
            res[:] = (disp.velocity - d)./sd
        end
    catch 
        map(swd.r) do res
            res[:] .= Inf
        end
    end
    chi2by2 = sum(reduce(vcat, swd.r).^2)
end

end