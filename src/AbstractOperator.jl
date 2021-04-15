module AbstractOperator
# basic operator types to do inversion with
abstract type Operator end
abstract type Operator1D <: Operator end
abstract type Operator2D <: Operator end

# any physics function must be a subtype of the
# above operators and needs to provide a misfit
function get_misfit end

# many geophysics data types require a sounding type
abstract type Sounding end
# these soundings often need to be fed into a function
# to make a physics operator with the relevant sounding data
function makeoperator end
# for nuisances this may change from sounding to sounding
function make_tdgp_opt end

export Operator, Operator1D, Operator2D
end
