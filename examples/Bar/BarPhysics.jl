# If you aren't trying to add to transD_GP but want to use what's
# already available.
module BarPhysics
# completely independent of transD_GP
# you can use whatever packages you want here

# this is probably where you'll store all variables needed for your physics
# e.g., wavenumbers and so on
mutable struct Bar
    x         :: Int
    somefield :: Vector{Float64}
end

# this is what will return the computation from our physics engine
function returnphysics!(F::Bar, m::AbstractArray)
    n = length(F.somefield)
    F.somefield .= F.x*rand(n) .+ sum(m[:]) # pre-allocated assignment
end

end