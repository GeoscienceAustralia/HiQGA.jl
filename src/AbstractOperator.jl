module AbstractOperator
# basic operator types to do inversion with
abstract type Operator end
abstract type Operator1D <: Operator end
abstract type Operator2D <: Operator end

# any physics function must be a subtype of the
# above operators and needs to provide a misfit
function get_misfit end

# for gradient based inversion
function getresidual end

# many geophysics data types require a sounding type
abstract type Sounding end
# we sometimes return things from a sounding in an inversion file
function returnforwrite end
# these soundings often need to be fed into a function
# to make a physics operator with the relevant sounding data
function makeoperator end
# for nuisances this may change from sounding to sounding
function make_tdgp_opt end
# to do AEM soundings in a loop
function loopacrossAEMsoundings end
# to plot forwards from a sounding
function plotmodelfield! end
# to get number of data from sounding
function getndata end
# to make bounds from sounding
function makebounds end
# to get optn with nuisances from existing struct
function getoptnfromexisting end

export Operator, Operator1D, Operator2D, Sounding
end
