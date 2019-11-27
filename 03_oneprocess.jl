# Reconstruct a badly sampled image with parsimonious representation
pwd()[2] == ':' && (ENV["MPLBACKEND"]="qt4agg")
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using PyPlot, Random, GP, TransD_GP, Revise,
        Distributed, Images, FileIO, DelimitedFiles
import MCMC_Driver, Plot2D
## get random points from image
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
## plot image
f1, ax1 = plt.subplots(1,2,figsize=(10,5), sharex=true, sharey=true)
im1 = ax1[1].imshow(f, extent=[x[1],x[end],y[end],y[1]])
cb1 = colorbar(im1, ax=ax1[1])
ax1[2].imshow(f, extent=[x[1],x[end],y[end],y[1]], alpha=0.0)
im2 = ax1[2].scatter(Xall[1,:], Xall[2,:], c=noisyd[:], s=10)
cb2 = colorbar(im2)
ax1[2].axis([x[1],x[end],y[end],y[1]])
MCMC_Driver.nicenup(gcf(), fsize=14)

## now for McMC inversion parameters
fdataname = "2Dtest_smooth"
nmin, nmax = 2, 20
λ, δ = [150.0, 150.0], 0.2
fbounds = [-1 2]
demean = true
sdev_prop = 0.1
sdev_pos = [10.0, 10.0]
pnorm = 2.
debug = false
MLnoise = true

xbounds = [x[1] x[end];y[1] y[end]]
## run McMC
opt_in = TransD_GP.Options(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = Xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        save_freq = 500,
                        dispstatstoscreen = false,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        debug = debug,
                        fdataname = fdataname
                        )

opt_EM_in  = MCMC_Driver.EMoptions(sd=δtry)
m = TransD_GP.init(opt_in)
m_true = deepcopy(m)
m_true.fstar[:] = f[:]
opt_EM_in.MLnoise = false
current_misfit = [MCMC_Driver.get_misfit(m_true, noisyd, opt_in, opt_EM_in)]
@info "True RMS error is" sqrt(2.0*current_misfit[1]/sum(.!(isnan.(noisyd))))

opt_EM_in.MLnoise = MLnoise
stat_in = TransD_GP.Stats()
wp_in = TransD_GP.open_history(opt_in)
Temp = 1.0
current_misfit = [MCMC_Driver.get_misfit(m, noisyd, opt_in, opt_EM_in)]
# run 4001 steps on single process - about 113 seconds
@info MCMC_Driver.get_misfit(m, noisyd, opt_in, opt_EM_in)
function foo()
    for isample=1:4001
        MCMC_Driver.do_mcmc_step(m, opt_in, stat_in, current_misfit, noisyd[:],
            Temp, isample, opt_EM_in, wp_in)
    end
    nothing
end
@time foo()
TransD_GP.close_history(wp_in)
@info MCMC_Driver.get_misfit(m, noisyd, opt_in, opt_EM_in)
figure()
imshow(reshape(m.fstar,length(y), length(x)), extent=[x[1],x[end],y[end],y[1]])