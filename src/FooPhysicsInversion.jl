# make sure to include this file in src/GeophysOperator.jl
# AFTER FooPhysics.jl has been included
module FooPhysicsInversion
using ..AbstractOperator
# using AbstractOperator (above) says that transD_GP can 
# use this kind of physics, through a misfit function

import ..AbstractOperator.get_misfit
# get_misfit (above) extends all the other misfit functions 
# in transD_GP. Because we are extending a method
# we must use import instead of using.

import ..Model, ..Options
# the above line now recognizes GP models and McMC options defined in src/TransD_GP_MCMC.jl

# The above three lines are the bare minimum we need to use transD_GP

# But wait! we need to specialize this inversion module 
# for FooPhysics ... so let's do that. This is basically the 
# geophysics "forward call"
using ..FooPhysics  
import ..FooPhysics.returnphysics!

# Now let's actually define our FooInversion operator and for simplicity
# say it is a subtype of 1D operator (operators are defined in src/AbstractOperator.jl)
mutable struct FooInversion<:Operator1D
    d # some observations data
    F # some forward physics struct
end   

# some constructor for FooInversion
function FooInversion() 
    FooInversion([0, 0.], FooPhysics.Foo(20, rand(1:100))) # return value
end

# Now let's define (extend) a misfit function for
# a 1D model parameter vector m, the physics for which 
# can be computed through methods in FooPhysics
# using an operator of type Foo, using McMC prior bounds 
# and options stored in a struct of type Options
function get_misfit(m::Model, opt::Options, foo::FooInversion)
    get_misfit(m.fstar, opt, foo) # dispatches to helper function below
end
# above defined function and type signature MUST be defined

function get_misfit(x::AbstractArray, opt::Options, foo::FooInversion)
    # m.fstar is the 1D array interpolated by the GP mean
    # at locations opt.xall
    returnphysics!(foo.F, x) # think of this as F(x)
    chi2by2 = 0*foo.d[1] - 0*foo.F.somefield[1]
end
# this misfit function returns the same value for very call 
# so we're sampling the prior
# this type signature is easier to debug than feeding in a GP model

end