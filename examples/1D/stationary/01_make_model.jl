# # Example - Bayesian nonlinear regression
# ## Setup
cd(@__DIR__)
using PyPlot, Random, LinearAlgebra, Distributed, DelimitedFiles, HiQGA.transD_GP
# Import data representing the "true" function
Random.seed!(200)
func_file = "../stationary/func2.txt"
x = readdlm(func_file, ',', Float64, '\n')[:,1]
y = readdlm(func_file, ',', Float64, '\n')[:,2]
σ, fractrain = 0.275, 0.5
ntrain = round(Int, fractrain*length(y))
ynoisy = similar(y) .+ NaN
linidx = randperm(length(y))[1:ntrain]
ynoisy[linidx] = y[linidx] + σ*randn(ntrain)
noisy_file = "../stationary/noisyjump1D.txt"
# done on Julia 1.8.1 to stop rng from changing -- no need to uncomment 
# as file with ynoisy is provided on GitHub, noisyjump1D.txt
# writedlm(noisy_file, ynoisy)
ynoisy = readdlm(noisy_file)[:]
line = transD_GP.Line(ynoisy;useML=false, σ=σ, calcjacobian=false)
figure(figsize=(4,3))
plot(x[:], y)
plot(x[:], ynoisy, ".m", alpha=0.5)
xlabel("x")
ylabel("f(x)")
plt.tight_layout()
gcf()
# Save the figure
savefig("jump1D.png", dpi=300)
