module HiQGA

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
