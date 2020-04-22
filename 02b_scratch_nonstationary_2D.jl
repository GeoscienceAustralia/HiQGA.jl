srcdir = "/Users/anray/Desktop/TDGP/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using PyPlot, Random, Statistics, Linear Algebra, Test, Revise,
MCMC_Driver, GP
##
Random.seed!(13)
xtest = collect(0:0.01:1)'
λ = [0.05]
δ = 0.02
δtry = 5δ
p = 2
demean = true

y = zeros(size(xtest,2))
change2 = round(Int, 3/4*length(y))
change1 = round(Int, 1/4*length(y))
y[change2:end] .= 2
y[change1:change2-1] .= 1
ynoisy = y + δ*randn(length(y))
dec = 4
ytrain = ynoisy[1:dec:end]
xtrain = xtest[1:dec:end]'

# f1 = figure()
# plot(xtest[:][:],y, "-k",label="true")
# plot(xtrain[:],ytrain, "+m", label="train",markersize=10)
#
ytest, σ2, σ_prior  = GP.GPfit(ytrain, xtrain, xtest, λ, δtry, p=p, demean=demean)
# plot(xtest[:],ytest,label="testing")
# plot(xtest[:],ytest+2sqrt.(diag(σ2)),"--r",alpha=0.2)
# plot(xtest[:],ytest-2sqrt.(diag(σ2)),"--b",alpha=0.2)

## non stationary
λ[1] = 0.05
λtest = λ[1]*ones(size(xtest))
λtrain = λ[1]*ones(Float64, size(xtrain))
ytest_ns, σ2_ns, σ_prior_ns  = GP.GPfit(ytrain, xtrain, xtest, λtest, λtrain, δtry, p=p, demean=demean)

f2, ax2 = plt.subplots(2,1, sharex=true, figsize=(7,7))
ax2[1].plot(xtest[:],y, "-k",label="true", linewidth=2)
ax2[1].plot(xtrain[:],ytrain, "xm", label="train",markersize=10, alpha=0.8)
ax2[1].plot(xtest[:],ytest_ns,label="fixed λ", linewidth=3)
# plot(xtest[:],ytest_ns+2sqrt.(diag(σ2_ns)),"--r",alpha=0.2)
# plot(xtest[:],ytest_ns-2sqrt.(diag(σ2_ns)),"--b",alpha=0.2)
s = gca()
##
lscale_max = 0.05
dropoff = 0.0025
λtest = lscale_max*(1 .+eps() .- (exp.(-1/dropoff*abs.(xtest.-xtest[change1+1]).^2) +
                    exp.(-1/dropoff*abs.(xtest.-xtest[change2-2]).^2)))
λtrain[:] = λtest[1:dec:end]
ytest_ns, σ2_ns, σ_prior_ns  = GP.GPfit(ytrain, xtrain,
    xtest, λtest, λtrain, δtry, p=2, demean=demean, nogetvars=false)

ax2[1].plot(xtest[:],ytest_ns,label="variable λ", linewidth=3, alpha=0.7, "-r")
ax2[2].plot(xtest[:],λtest[:], linewidth=3, color="r", alpha=0.7)
ax2[1].set_ylabel(L"f(y)")
ax2[2].set_xlabel(L"y")
ax2[2].set_ylabel(L"\lambda")
ax2[1].grid()
ax2[2].grid()
MCMC_Driver.nicenup(gcf(), fsize=14)
savefig("variable.png", dpi=300)
## timing
@info "meshkernel with broadcast"
@time for i = 1:1000
    GP.meshkernel(xtrain, xtest, λtest, λtrain, p)
end
@info "mapkernel with pairwise"
@time for i = 1:1000
    GP.mapkernel(xtrain, xtest, λtest, λtrain)
end
@info "loops and full expression"
@time for i = 1:1000
    GP.makekernel(xtrain, xtest, λtest, λtrain, p)
end
## unit tests
@testset "Non stationary kernel tests" begin
    @testset "broadcast and map" begin
        @test norm(mean(GP.meshkernel(xtrain, xtest, λtest, λtrain, p) -
                    GP.mapkernel(xtrain, xtest, λtest, λtrain))) < 1e-12
        end
    @testset "map and full kernel" begin
        @test norm(mean(GP.makekernel(xtrain, xtest, λtest, λtrain, p) -
                GP.mapkernel(xtrain, xtest, λtest, λtrain)')) < 1e-12
        end
end
