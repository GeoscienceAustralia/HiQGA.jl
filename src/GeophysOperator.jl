import .AbstractOperator.Operator
import .AbstractOperator.get_misfit
import .AbstractOperator.Sounding
import .AbstractOperator.makeoperator
import .AbstractOperator.make_tdgp_opt
include("CommonToAll.jl")
using .CommonToAll
include("LineRegression.jl")
using .LineRegression
include("SurfaceRegression.jl")
using .SurfaceRegression
include("CSEM1DEr.jl") # add this before using HED for CSEM
using .CSEM1DEr
include("CSEM1DInversion.jl")
using .CSEM1DInversion
include("ImageRegression.jl")
using .ImageRegression
include("AEM_VMD_HMD.jl") # add this before using VMD or HMD for land/AEM geophysics
using .AEM_VMD_HMD
include("SkyTEM1DInversion.jl")
using .SkyTEM1DInversion
include("TEMPEST1DInversion.jl")
using .TEMPEST1DInversion
export Operator, get_misfit, Sounding, makeoperator, make_tdgp_opt
