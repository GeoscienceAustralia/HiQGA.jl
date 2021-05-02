using PyPlot, Random, Revise, Statistics, LinearAlgebra,
      Distributed, DelimitedFiles, transD_GP
## 1D functions
Random.seed!(200)
x = readdlm("func2.txt", ',', Float64, '\n')[:,1]
y = readdlm("func2.txt", ',', Float64, '\n')[:,2]
σ, fractrain = 0.275, 0.5
ntrain = round(Int, fractrain*length(y))
ynoisy = similar(y) .+ NaN
linidx = randperm(length(y))[1:ntrain]
ynoisy[linidx] = y[linidx] + σ*randn(ntrain)
line = transD_GP.Line(ynoisy;useML=false, σ=σ)
figure(figsize=(4,3))
plot(x[:], y)
plot(x[:], ynoisy, ".m", alpha=0.5)
xlabel("x")
ylabel("f(x)")
plt.tight_layout()
savefig("jump1D.png", dpi=300)
