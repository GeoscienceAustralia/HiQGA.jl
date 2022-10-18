module HiQGA
module transD_GP
    include("GP.jl")
    include("TransD_GP_MCMC.jl")
    include("AbstractOperator.jl")
    include("MCMC_Driver.jl")
    include("GradientInversion.jl")
    include("GeophysOperator.jl")
    include("AEMnoNuisanceGradientInversionTools.jl")
    using .AEMnoNuisanceGradientInversionTools
end
export transD_GP
end