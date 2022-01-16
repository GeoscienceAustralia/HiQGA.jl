using transD_GP
##
ntimesperdecade = 2
T = 10 .^(-3:1/ntimesperdecade:6)
f = 1 ./T
## model 1 conductor above
ρ = [10, 100]
z = [0, 5000]
## make a synthetic
F = transD_GP.MT1DInversion.create_synthetic(ρ, z, f, rseed=1)
transD_GP.get_misfit(F.d_log10_ρ, F.d_phase_deg, F.σ_log10_ρ, F.σ_phase_deg, f, ρ, z)