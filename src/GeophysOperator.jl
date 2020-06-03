module GeophysOperator
import AbstractOperator.Operator
import AbstractOperator.get_misfit
include("LineRegression.jl")
using .LineRegression
include("CommonToAll.jl")
using .CommonToAll
#include("ImageRegression.jl")
export Operator, get_misfit
end
