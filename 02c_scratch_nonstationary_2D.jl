srcdir = "/Users/anray/Desktop/TDGP/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using PyPlot, Random, Statistics, LinearAlgebra, Test, Revise,
MCMC_Driver, GP
##
Random.seed!(11)
xtest = collect(0:0.005:1)'
λ² = [0.1].^2
δ = 0.01
δtry = 10δ
p = 2
demean = true
K = GP.Mat32()

f = zeros(size(xtest,2))
change1 = round(Int, 1/2*length(f))
f[change1:end] .= 1
f .-= 0.5
ynoisy = f + δ*randn(length(f))
dec = 8
ytrain = ynoisy[1:dec:end]'
xtrain = xtest[1:dec:end]'

# f1 = figure()
# plot(xtest[:][:],f, "-k",label="true")
# plot(xtrain[:],ytrain, "+m", label="train",markersize=10)
#
ytest, σ2, σ_prior  = GP.GPfit(K, ytrain, xtrain, xtest, λ², δtry, p=p, demean=demean)
# plot(xtest[:],ytest,label="testing")
# plot(xtest[:],ytest+2sqrt.(diag(σ2)),"--r",alpha=0.2)
# plot(xtest[:],ytest-2sqrt.(diag(σ2)),"--b",alpha=0.2)

## non stationary
λ²test = λ²[1]*ones(size(xtest))
λ²train = λ²[1]*ones(Float64, size(xtrain))
ytest_ns, σ2_ns, σ_prior_ns  = GP.GPfit(K, ytrain, xtrain, xtest, λ²test, λ²train, δtry, p=p, demean=demean)
@test (norm(mean(ytest - ytest_ns)) < 1e-12)
f2, ax2 = plt.subplots(2,1, sharex=true, figsize=(7,7))
ax2[1].plot(xtest[:],f, "-k",label="true", linewidth=2)
ax2[1].plot(xtrain[:],ytrain, "xm", label="train",markersize=10, alpha=0.8)
ax2[1].plot(xtest[:],ytest_ns,label="fixed λ", linewidth=3)
# plot(xtest[:],ytest_ns+2sqrt.(diag(σ2_ns)),"--r",alpha=0.2)
# plot(xtest[:],ytest_ns-2sqrt.(diag(σ2_ns)),"--b",alpha=0.2)
##
lscale_max² = λ²[1]
dropoff = 0.00025
λ²test = lscale_max²*(1 .+eps() .- exp.(-1/dropoff*abs.(xtest.-xtest[change1+1]).^2)).^2
λ²train[:] = λ²test[1:dec:end]
ytest_ns, σ2_ns, σ_prior_ns  = GP.GPfit(K, ytrain, xtrain,
    xtest, λ²test, λ²train, δtry, p=2, demean=demean, nogetvars=false)

ax2[1].plot(xtest[:],ytest_ns,label="variable λ", linewidth=3, alpha=0.7, "-r")
ax2[2].plot(xtest[:],sqrt.(λ²test[:]), linewidth=3, color="r", alpha=0.7)
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
    GP.meshkernel(K, xtrain, xtest, λ²test, λ²train, p)
end
@info "mapkernel with pairwise"
@time for i = 1:1000
    GP.mapkernel(K, xtrain, xtest, λ²test, λ²train)
end
@info "loops and full expression"
@time for i = 1:1000
    GP.makekernel(K, xtrain, xtest, λ²test, λ²train, p)
end
## unit tests
@testset "Non stationary kernel tests" begin
    @testset "broadcast and map" begin
        @test norm(mean(GP.meshkernel(K, xtrain, xtest, λ²test, λ²train, p) -
                    GP.mapkernel(K, xtrain, xtest, λ²test, λ²train))) < 1e-12
        end
    @testset "map and full kernel" begin
        @test norm(mean(GP.makekernel(K, xtrain, xtest, λ²test, λ²train, p) -
                GP.mapkernel(K, xtrain, xtest, λ²test, λ²train)')) < 1e-12
        end
end
## now to 2D
Random.seed!(11)
x, y = xtest, xtest[1:round(Int, length(xtest)/2)]
x2D_test = zeros(2,length(x)*length(y))
λ²2D_test = zeros(2,length(x)*length(y))
for i in 1:size(x2D_test,2)
    xid, yid = Tuple(CartesianIndices((length(x),length(y)))[i])
    x2D_test[:,i] = [x[xid]; y[yid]]
    λ²xid, λ²yid = Tuple(CartesianIndices((length(x),length(y)))[i])
    λ²2D_test[1,i] = λ²test[λ²xid]
end
λ²2D_test[2,:] .= abs(y[end]-y[1])^2
y2D       = [f[i] for i in 1:length(x), j in 1:length(y)]
y2D_noisy = y2D + δ*randn(size(y2D))
y2D_train = y2D_noisy[1:dec:end, 1:dec:end][:]'
lidx = LinearIndices((1:length(x),1:length(y)))[1:dec:end, 1:dec:end][:]
λ²2D_train = λ²2D_test[:,lidx]
x2D_train = x2D_test[:,lidx]
## fit GPs
y2D_test_ns, = GP.GPfit(K, y2D_train, x2D_train,
    x2D_test, λ²2D_test, λ²2D_train, δtry, p=2, demean=demean, nogetvars=true)
y2D_test, = GP.GPfit(K, y2D_train, x2D_train,x2D_test,
            [λ²[1], λ²2D_test[2,1]],
            δtry, p=2, demean=demean, nogetvars=true)
## plot
fig, ax = plt.subplots(1,3, sharex=true, sharey=true)
im1 = ax[1].imshow(y2D_noisy)
ax[1].set_title("data");
colorbar(im1, ax=ax[1])
im2 = ax[2].imshow(reshape(y2D_test, size(y2D, 1), size(y2D, 2)))
ax[2].set_title("stationary GP")
colorbar(im2, ax=ax[2])
im3 = ax[3].imshow(reshape(y2D_test_ns, size(y2D, 1), size(y2D, 2)))
ax[3].set_title("non stationary GP")
colorbar(im3, ax=ax[3])
## see how fast you can change a row
Kstar = GP.mapkernel(K, x2D_train, x2D_test, λ²2D_test, λ²2D_train)
K_y  = GP.mapkernel(K, x2D_train, x2D_train, λ²2D_train, λ²2D_train) +
    Matrix(δtry^2*I, (size(x2D_train,2)), (size(x2D_train,2)))
A = copy(Kstar)
B = copy(K_y)
# change a training point
n = rand(1:length(y2D_train))
x2D_train[:,n] = x2D_test[:,n]
λ²2D_train[:,n] = λ²2D_test[:,n]
y2D_train[n]   = y2D_noisy[CartesianIndices((length(x),length(y)))[n]]
D, = GP.GPfit(K, y2D_train, x2D_train,
    x2D_test, λ²2D_test, λ²2D_train, δtry, p=2, demean=demean, nogetvars=true)
@time begin
Kstarv = @view Kstar[:,n]
map!(x²->x²,Kstarv,GP.colwise(K, x2D_train[:,n], x2D_test, λ²2D_train[:,n], λ²2D_test))
K_yv = @view K_y[n,:]
map!(x²->x²,K_yv,GP.colwise(K, x2D_train[:,n], x2D_train, λ²2D_train[:,n], λ²2D_train))
K_y[:,n] = K_y[n,:]
K_y[n,n] = K_y[n,n] + δtry^2
my = mean(y2D_train, dims=2)
y2D_train = y2D_train .- my
U = cholesky(K_y).U
C = my' .+ Kstar*(U\(U'\y2D_train'))
end
@test norm(mean(C - D)) < 1e-12
