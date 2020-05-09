using Distances, NearestNeighbors
mns = TransD_GP.init(opt, log10λ)
mnsold = deepcopy(mns)
log10λold = deepcopy(log10λ)
TransD_GP.birth!(log10λ, optlog10λ, mns, opt.δ)
TransD_GP.sync_model!(mns, opt)
λ² = log10λ.fstar
idxs = TransD_GP.gettrainidx(opt.kdtree, mnsold.xtrain, mnsold.n)
ftest, = GP.GPfit(K, mnsold.ftrain[:,1:mnsold.n], mnsold.xtrain[:,1:mnsold.n],
        opt.xall, λ², λ²[:,idxs], opt.δ, p=2, demean=demean, nogetvars=true)
##
m = deepcopy(log10λold)
m.fstar = 0.5log10.(m.fstar)
fig,ax = plt.subplots(1,2, sharex=true, sharey=true)
vmin, vmax = minimum(m.fstar[1,:]), maximum(m.fstar[1,:])
im1 = ax[1].imshow(reshape(m.fstar[1,:],length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
    vmin=vmin, vmax=vmax)
colorbar(im1, ax=ax[1])
ax[1].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],marker="+",c="r")
ax[1].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],c=m.ftrain[1,1:m.n], alpha=0.8,
vmin=vmin, vmax=vmax)
vmin, vmax = minimum(m.fstar[2,:]), maximum(m.fstar[2,:])
im2 = ax[2].imshow(reshape(m.fstar[2,:],length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
    vmin=vmin, vmax=vmax)
colorbar(im2, ax=ax[2])
ax[2].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],marker="+",c="r")
ax[2].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],c=m.ftrain[2,1:m.n], alpha=0.8,
    vmin=vmin, vmax=vmax)
##
fig,ax = plt.subplots(1,2, sharex=true, sharey=true)
vmin, vmax = minimum(mold.fstar), maximum(mold.fstar)
im1 = ax[1].imshow(reshape(mold.fstar,length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
    vmin=vmin, vmax=vmax)
colorbar(im1, ax=ax[1])
ax[1].scatter(mold.xtrain[1,1:mold.n], mold.xtrain[2,1:mold.n],marker="+",c="r")
ax[1].scatter(mold.xtrain[1,1:mold.n], mold.xtrain[2,1:mold.n],c=mold.ftrain[1:mold.n], alpha=0.8)
vmin, vmax = minimum(mns.fstar), maximum(mns.fstar)
im2 = ax[2].imshow(reshape(mns.fstar,length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
    vmin=vmin, vmax=vmax)
colorbar(im2, ax=ax[2])
ax[2].scatter(mns.xtrain[1,1:mns.n], mns.xtrain[2,1:mns.n],marker="+",c="r")
ax[2].scatter(mns.xtrain[1,1:mns.n], mns.xtrain[2,1:mns.n],c=mns.ftrain[1:mns.n], alpha=0.8,
    vmin=vmin, vmax=vmax)
