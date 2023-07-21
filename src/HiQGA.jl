module HiQGA

# using PyCall 
# using PyPlot

# const matplotlib = PyNULL()
# const plt = PyNULL()

# function __init__()
#     @info "I made it here"
#     copy!(matplotlib, pyimport_conda("matplotlib", "matplotlib", "conda-forge"))
#     copy!(plt, pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge"))
# end

module transD_GP
    include("GP.jl")
    include("TransD_GP_MCMC.jl")
    include("AbstractOperator.jl")
    include("MCMC_Driver.jl")
    include("GradientInversion.jl")
    include("SoundingDistributor.jl")
    using .SoundingDistributor
    include("GeophysOperator.jl")
    include("AEMnoNuisanceGradientInversionTools.jl")
    using .AEMnoNuisanceGradientInversionTools
    include("AEMwithNuisanceGradientInversionTools.jl")
    using .AEMwithNuisanceGradientInversionTools
    include("AEMnoNuisanceMcMCInversionTools.jl")
    using .AEMnoNuisanceMcMCInversionTools
    include("AEMwithNuisanceMcMCInversionTools.jl")
    using .AEMwithNuisanceMcMCInversionTools
end
export transD_GP
end
