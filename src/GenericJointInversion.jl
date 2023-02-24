module GenericJointInversion

using ..AbstractOperator
import ..AbstractOperator.get_misfit
import ..Model, ..Options

struct JointOperator <: Operator
    Fvec :: Vector{x} where x<:Operator
end    

function get_misfit(m::Model, opt::Options, F::JointOperator)
    chi2by2 = 0.0
    for Fop in F.Fvec
        chi2by2 += get_misfit(m, opt, Fop) # dispatches to specific F
    end
    chi2by2    
end

end