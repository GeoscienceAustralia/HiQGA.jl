module MT1DInversion
using ..AbstractOperator
import ..AbstractOperator.get_misfit
import ..Model, ..Options
using ..MT1D
using Random

mutable struct MT1DZ <: Operator1D
    d_log10_ρ
    d_phase_deg
    σ_log10_ρ
    σ_phase_deg
    freqs
end    

function get_misfit(m::Model, opt::Options, F::MT1DZ)
   
end
# above defined function and type signature MUST be defined

function get_misfit(ρ::AbstractArray, opt::Options, F::MT1DZ)
    
end

function get_misfit(d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg, freqs, ρ, z)
    Z = get_Z(freqs, ρ, z)
    ρₐ = MT1D.ρapp(freqs, Z)
    ϕ  = MT1D.phase(Z)
    r = (d_log10_ρ - log10.(ρₐ))./σ_log10_ρ + (d_phase_deg - ϕ)./σ_phase_deg
    r'r/2
end

get_Z(freqs, ρ, z) = MT1D.Z_f(freqs, ρ, diff(z))

function add_noise(ρ, z, freqs; noisefrac=0.05, rseed=1)
    Random.seed!(rseed)
    Z = get_Z(freqs, ρ, z)
    @show σ_real = noisefrac*abs.(Z)/2
    noise = σ_real.*( randn(length(freqs)) + 1im*randn(length(freqs)) )
    Z + noise
end

function create_synthetic(ρ, z, freqs; noisefrac=0.05, rseed=1)
    Znoisy = add_noise(ρ, z, freqs, noisefrac=noisefrac, rseed=rseed)
    d_log10_ρ = log10.(MT1D.ρapp(freqs, Znoisy))   
    d_phase_deg = MT1D.phase(Znoisy)
    σ_log10_ρ = noisefrac/log(10)
    σ_phase_deg = rad2deg(noisefrac/2)
    MT1DZ(d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg, freqs)
end

end
