module GeophysOperator
import AbstractOperator.Operator
import AbstractOperator.get_misfit
include("LineRegression.jl")
using .LineRegression
include("CommonToAll.jl")
using .CommonToAll
include("CSEM1DInversion.jl")
using .CSEM1DInversion
include("ImageRegression.jl")
using .ImageRegression
include("SkyTEM1DInversion.jl")
using .SkyTEM1DInversion
export Operator, get_misfit
end
