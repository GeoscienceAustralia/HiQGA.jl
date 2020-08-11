using PyPlot, Random, Revise, Statistics, LinearAlgebra,
      Distributed, DelimitedFiles
srcdir = dirname(dirname(dirname(pwd())))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using GP, TransD_GP, GeophysOperator, MCMC_Driver
##1D functions
Random.seed!(3)
x, y = GeophysOperator.linetestfunction()
σ = 0.05
ynoisy = σ*randn(size(y)) + y
line = GeophysOperator.Line(ynoisy;useML=false, σ=σ)
figure(figsize=(4,3))
plot(x[:], y)
plot(x[:], ynoisy, ".m", alpha=0.5)
savefig("jump1D.png", dpi=300)
