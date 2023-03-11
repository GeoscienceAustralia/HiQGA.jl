## get nonstat stuff first
# plot a mean and sdev
close("all")
using Statistics
M = transD_GP.assembleTat1(opt, :fstar, temperaturenum=1)
mns = reshape(mean(M), length(img.y), length(img.x))
stdns = reshape(std(M), length(img.y), length(img.x))
f, ax = plt.subplots(1,2, sharex=true, sharey=true, figsize=(9,4))
im1 = ax[1].imshow(mns, extent=[img.x[1],img.x[end],img.y[end],img.y[1]])
im2 = ax[2].imshow(stdns, extent=[img.x[1],img.x[end],img.y[end],img.y[1]], cmap="bone_r")
map(x->colorbar(x), (im1, im2))

## sections through the image
rownum, colnum = 195, 85
transD_GP.plot_image_posterior(opt, optlog10λ, img, burninfrac=0.5, rownum=rownum, colnum=colnum, nbins=100)

## plot stats, pixel values
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
linidx = .!isnan.(img.d)
r = vec(img.d[linidx] - img.f[linidx])
χ² = length(r)
ax = gcf().axes
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
## stats on length scales
transD_GP.getchi2forall(optlog10λ, fsize=8, alpha=0.5)
ax = gcf().axes
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")

