using PyPlot, Random, Revise, Statistics, LinearAlgebra,
      Distributed, DelimitedFiles
srcdir = dirname(dirname(dirname(pwd())))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using GP, TransD_GP, GeophysOperator, MCMC_Driver
##1D functions
Random.seed!(200)
x = readdlm("func2.txt", ',', Float64, '\n')[:,1]
y = readdlm("func2.txt", ',', Float64, '\n')[:,2]
σ = 0.55
dec = 1
till = round(Int, length(x)/2)
ynoisy = σ*randn(size(y)) + y
keep = copy(ynoisy)
ynoisy[till+1:end] .= NaN
ynoisy[till+1:dec:end] .= keep[till+1:dec:end]
line = GeophysOperator.Line(ynoisy;useML=false, σ=σ)
figure(figsize=(4,3))
plot(x[:], y)
plot(x[:], ynoisy, ".m", alpha=0.5)
xlabel("x")
ylabel("f(x)")
plt.tight_layout()
savefig("jump1D.png", dpi=300)
