srcdir = "/Users/anray/Desktop/TDGP/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using PyPlot, Random, Statistics, LinearAlgebra, Test, Revise,
MCMC_Driver, GP
##
Random.seed!(13)
xtest = collect(0:0.005:1)'
λ = [0.05]
δ = 0.02
δtry = 5δ
p = 2
demean = true

f = zeros(size(xtest,2))
change2 = round(Int, 3/4*length(f))
change1 = round(Int, 1/4*length(f))
f[change2:end] .= 2
f[change1:change2-1] .= 1
ynoisy = f + δ*randn(length(f))
dec = 4
ytrain = ynoisy[1:dec:end]
xtrain = xtest[1:dec:end]'

# f1 = figure()
# plot(xtest[:][:],f, "-k",label="true")
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
ax2[1].plot(xtest[:],f, "-k",label="true", linewidth=2)
ax2[1].plot(xtrain[:],ytrain, "xm", label="train",markersize=10, alpha=0.8)
ax2[1].plot(xtest[:],ytest_ns,label="fixed λ", linewidth=3)
# plot(xtest[:],ytest_ns+2sqrt.(diag(σ2_ns)),"--r",alpha=0.2)
# plot(xtest[:],ytest_ns-2sqrt.(diag(σ2_ns)),"--b",alpha=0.2)
##
lscale_max = 0.05
dropoff = 0.0005
λtest = lscale_max*(1 .+eps() .- (exp.(-1/dropoff*abs.(xtest.-xtest[change1+1]).^2) +
                    exp.(-1/dropoff*abs.(xtest.-xtest[change2-2]).^2)))
λtrain[:] = λtest[1:dec:end]
ytest_ns, σ2_ns, σ_prior_ns  = GP.GPfit(ytrain, xtrain,
    xtest, λtest, λtrain, δtry, p=2, demean=demean, nogetvars=false)

ax2[1].plot(xtest[:],ytest_ns,label="variable λ", linewidth=3, alpha=0.7, "-r")
ax2[2].plot(xtest[:],λtest[:], linewidth=3, color="r", alpha=0.7)
ax2[1].set_ylabel(L"f(x)")
ax2[2].set_xlabel(L"x")
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
## now to 2D
x, y = xtest, xtest[1:round(Int, length(xtest)/2)]
x2D_test = zeros(2,length(x)*length(y))
λ2D_test = zeros(2,length(x)*length(y))
for i in 1:size(x2D_test,2)
    xid, yid = Tuple(CartesianIndices((length(x),length(y)))[i])
    x2D_test[:,i] = [x[xid]; y[yid]]
    λxid, λyid = Tuple(CartesianIndices((length(x),length(y)))[i])
    λ2D_test[1,i] = λtest[λxid]
end
λ2D_test[2,:] .= abs(y[end]-y[1])
y2D       = [f[i] for i in 1:length(x), j in 1:length(y)]
y2D_noisy = y2D + δ*randn(size(y2D))
y2D_train = y2D_noisy[1:dec:end, 1:dec:end][:]
lidx = LinearIndices((1:length(x),1:length(y)))[1:4:end, 1:4:end][:]
λ2D_train = λ2D_test[:,lidx]
x2D_train = x2D_test[:,lidx]
y2D_test_ns, = GP.GPfit(y2D_train, x2D_train,
    x2D_test, λ2D_test, λ2D_train, δtry, p=2, demean=demean, nogetvars=true)
y2D_test, = GP.GPfit(y2D_train, x2D_train,
        x2D_test, [λ[1], abs(y[end]-y[1])], δtry, p=2, demean=demean, nogetvars=true)
