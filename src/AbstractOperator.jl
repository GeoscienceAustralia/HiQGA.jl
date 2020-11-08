module AbstractOperator
abstract type Operator end
abstract type Operator1D <: Operator end
abstract type Operator2D <: Operator end
function get_misfit end
export Operator, Operator1D, Operator2D
end
