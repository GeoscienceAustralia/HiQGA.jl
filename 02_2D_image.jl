pwd()[2] == ':' && (ENV["MPLBACKEND"]="qt4agg")
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using PyPlot, Random, GP, JLD, TransD_GP, Revise, Distributed, Images, FileIO
import MCMC_Driver, Plot2D
##
Random.seed!(12)
sd = 0.05
fractrain = 0.02
dec = 2
f = Gray.(load("4.2.01.tiff"))
f = convert(Array{Float64, 2}, f)[1:dec:end,:1:dec:end]
f = -1 .+ 3*f
f = imfilter(f,Kernel.gaussian(7))
dx = 10.0
x = 0:dx:dx*size(f,2)-1
y = 0:dx:dx*size(f,1)-1
λ = [150.0; 150.0]
δtry = sd*max(f...)
noisyd = NaN .+ zeros(Float64, size(f))
ntrain = round(Int, fractrain*length(f))
##
Xtrain = zeros(2,0)
linidx = randperm(length(f))[1:ntrain]
lgood = zeros(Int, 0)
for (i,l) in enumerate(linidx)
    row, col = Tuple(CartesianIndices(f)[l])
    if y[row] > 1000
        noisyd[row, col] = f[row, col] + δtry*randn()
        push!(lgood, l)
        global Xtrain = hcat(Xtrain, [x[col]; y[row]])
    else
        if rem(row,4) == 0 && rem(col,4) == 0
            noisyd[row, col] = f[row, col] + δtry*randn()
            push!(lgood, l)
            global Xtrain = hcat(Xtrain, [x[col]; y[row]])
        end
    end
end
ftrain = noisyd[lgood]
Xall = zeros(2,length(noisyd))
for i in 1:length(noisyd)
    yid, xid = Tuple(CartesianIndices(f)[i])
    Xall[:,i] = [x[xid]; y[yid]]
end
##
#f1, ax1 = plt[:subplots](1,2,figsize=(10,5), sharex=true, sharey=true)
#im1 = ax1[1][:imshow](f, extent=[x[1],x[end],y[end],y[1]])
#cb1 = colorbar(im1, ax=ax1[1])
#ax1[2][:imshow](f, extent=[x[1],x[end],y[end],y[1]], alpha=0.0)
#im2 = ax1[2][:scatter](Xall[1,:], Xall[2,:], c=noisyd[:], s=10)
#cb2 = colorbar(im2)
#ax1[2][:axis]([x[1],x[end],y[end],y[1]])
## ax1[2][:invert_yaxis]()
#MCMC_Driver.nicenup(gcf(), fsize=14)
#savefig("2D_setup.png", dpi=300)
#@time ftest, = GP.GPfit(ftrain, Xtrain, Xall, λ, δtry, nogetvars=true)
#figure()
#imshow(reshape(ftest,length(y), length(x)), extent=[x[1],x[end],y[end],y[1]])
#colorbar()
## scatter(Xtrain[1,:], Xtrain[2,:], c=ftrain, s=50)
## now for inversion parameters
fdataname = "2Dtest_smooth"
nmin, nmax = 2, 100
λ, δ = [150.0, 150.0], 0.2
fbounds = [-1 2]
demean = true
sdev_prop = 0.1
sdev_pos = [10.0, 10.0]
pnorm = 2.
debug = false
MLnoise = true

xbounds = [x[1] x[end];y[1] y[end]]
opt_in = TransD_GP.Options(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = Xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        save_freq = 100,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        debug = debug,
                        costs_filename = "misfits_"*fdataname*".bin",
                        fstar_filename = "models_"*fdataname*".bin",
                        x_ftrain_filename = "points_"*fdataname*".bin"
                        )

opt_EM_in  = MCMC_Driver.EMoptions(sd=δtry)
m_true = TransD_GP.init(opt_in)
m_true.fstar[:] = f[:]
opt_EM_in.MLnoise = false
@info "RMS error is" sqrt(2.0*MCMC_Driver.get_misfit(m_true, noisyd, opt_in, opt_EM_in)/sum(.!(isnan.(noisyd))))
opt_EM_in.MLnoise = MLnoise
## run
nsamples = 10001
nchains = 2
Tmax = 2.5
rmprocs(workers()); addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed, Revise
@everywhere import MCMC_Driver
# m, opt, stat, opt_EM, d, current_misfit = MCMC_Driver.init_chain_darrays(opt_in, opt_EM_in, noisyd[:])
##
@time begin
    misfit, T0loc = MCMC_Driver.main(opt_in, noisyd, Tmax, nsamples, opt_EM_in)
end
save("misfit_T0_"*fdataname*".jld", "misfit", misfit, "T0loc", T0loc)
##
burnin = 5000
M = TransD_GP.history(opt_in, stat=:fstar)
n = TransD_GP.history(opt_in, stat=:nodes)
x_ft = TransD_GP.history(opt_in, stat=:x_ftrain)
iter = TransD_GP.history(opt_in, stat=:iter)
misfit = load("misfit_T0_"*fdataname*".jld", "misfit")
T0loc = load("misfit_T0_"*fdataname*".jld", "T0loc")
f2, ax2 = plt[:subplots](2,1, sharex=true, figsize=(8,4))
ax2[1][:plot](iter,n)
ax2[1][:grid]()
ax2[2][:plot](iter, misfit)
ax2[2][:grid]()
ax2[2][:set_xlabel]("iterations")
ax2[2][:set_ylabel]("-log likelihood")
ax2[1][:set_ylabel]("# training")
MCMC_Driver.nicenup(gcf(), fsize=14)
# gca()[:get_legend]()[:remove]()
savefig("2D_conv.png", dpi=300)
#gca()[:set_ylim](50, 200)
# subplot(313)
# plot(iter, T0loc)
s = zeros(size(M[1]))
iburn = findfirst(iter.>burnin)
for i = iburn:length(M)
    global s+= M[i]
end
s = s/(1-iburn+length(M))
m = deepcopy(m_true)
m.fstar[:] = M[end]
opt_EM_in.MLnoise = false
@info "RMS error is" sqrt(2.0*MCMC_Driver.get_misfit(m, noisyd, opt_in, opt_EM_in)/sum(.!(isnan.(noisyd))))
opt_EM_in.MLnoise = MLnoise
f3, ax3 = plt[:subplots](1,2,figsize=(10,5), sharex=true, sharey=true)
nmodel = iburn
im1 = ax3[1][:imshow](reshape(M[nmodel],length(y), length(x)), extent=[x[1],x[end],y[end],y[1]])
ax3[1][:scatter](x_ft[nmodel][1:n[nmodel],1], x_ft[nmodel][1:n[nmodel],2], s=20, color="black", alpha=0.5)
cb1 = colorbar(im1, ax=ax3[1])
im2 = ax3[2][:imshow](reshape(s,length(y), length(x)), extent=[x[1],x[end],y[end],y[1]])
cb2 = colorbar(im2, ax=ax3[2])
MCMC_Driver.nicenup(gcf(), fsize=14)
# gca()[:get_legend]()[:remove]()
savefig("2D_final.png", dpi=300)
##
