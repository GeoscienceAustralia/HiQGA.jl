module GeophysOperator
import AbstractOperator.Operator
import AbstractOperator.get_misfit
include("LineRegression.jl")
using .LineRegression
include("CommonToAll.jl")
using .CommonToAll
#include("ImageRegression.jl")
export Operator, Operator1D, Operator2D, get_misfit
end
