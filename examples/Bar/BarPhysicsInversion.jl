# make sure to include this file 
# AFTER BarPhysics.jl has been included
module BarPhysicsInversion
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
using ..BarPhysics  
import ..BarPhysics.returnphysics!

# Now let's actually define our BarInversion operator and for simplicity
# say it is a subtype of 1D operator (operators are defined in src/AbstractOperator.jl)
mutable struct BarInversion<:Operator1D
    d # some observations data
    F # some forward physics struct
end   

# some constructor for BarInversion
function BarInversion() 
    BarInversion([0, 0.], BarPhysics.Bar(20, rand(1:100))) # return value
end

# Now let's define (extend) a misfit function for
# a 1D model parameter vector m, the physics for which 
# can be computed through methods in BarPhysics
# using an operator of type Bar, using McMC prior bounds 
# and options stored in a struct of type Options
function get_misfit(m::Model, opt::Options, Bar::BarInversion)
    get_misfit(m.fstar, opt, Bar) # dispatches to helper function below
end
# above defined function and type signature MUST be defined

function get_misfit(x::AbstractArray, opt::Options, Bar::BarInversion)
    # m.fstar is the 1D array interpolated by the GP mean
    # at locations opt.xall
    returnphysics!(Bar.F, x) # think of this as F(x)
    chi2by2 = 0*Bar.d[1] - 0*Bar.F.somefield[1]
end
# this misfit function returns the same value for very call 
# so we're sampling the prior
# this type signature is easier to debug than feeding in a GP model

end