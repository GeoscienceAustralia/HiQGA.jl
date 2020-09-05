## get nonstat stuff first
cd("../nonstat")
include("01_make_model.jl")
include("02_make_options.jl")
close("all")
using Statistics
M = GeophysOperator.assembleTat1(opt, :fstar, temperaturenum=1)
mns = reshape(mean(M), length(img.y), length(img.x))
stdns = reshape(std(M), length(img.y), length(img.x))
##
rownum, colnum = 195, 85
GeophysOperator.plot_image_posterior(opt, optlog10λ, img, burninfrac=0.5, rownum=rownum, colnum=colnum, nbins=100)
figure(1)
ax = gcf().axes
ax[7].plot(img.x, img.f[rownum,:], "--y", alpha=0.9)
savefig("post_rows_ns.png", dpi=300)
figure(2)
ax = gcf().axes
ax[7].plot(img.f[:,colnum], img.y,"--y", alpha=0.9)
savefig("post_col_ns.png", dpi=300)
# sampling stats
## plot stats
GeophysOperator.getchi2forall(opt, fsize=8, alpha=0.5)
gcf().text(0.02, 0.9, "a.", fontsize=14, color="red")
ax = gcf().axes;
linidx = .!isnan.(img.d)
r = vec(img.d[linidx] - img.f[linidx])
χ² = length(r)#r'*r/img.σ^2
figure(1)
ax[2].set_ylim(χ²/2 - 100, χ²/2 + 300)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
figure(2)
ax = gcf().axes
ax[2].set_ylim(χ²/2 - 100, χ²/2 + 300)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
savefig("img_conv_ns_1.png", dpi=300)
close("all")
GeophysOperator.getchi2forall(optlog10λ, fsize=8, alpha=0.5)
gcf().text(0.02, 0.9, "b.", fontsize=14, color="red")
figure(1)
ax = gcf().axes
ax[2].set_ylim(χ²/2 - 100, χ²/2 + 300)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
figure(2)
ax = gcf().axes
ax[2].set_ylim(χ²/2 - 100, χ²/2 + 300)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
savefig("img_conv_ns_2.png", dpi=300)
close("all")
## change to stat directory
cd("../stat")
include("../stat/02_make_options.jl")
close("all")
M = GeophysOperator.assembleTat1(opt, :fstar, temperaturenum=1)
m = reshape(mean(M), length(img.y), length(img.x))
stds = reshape(std(M), length(img.y), length(img.x))
GeophysOperator.plot_image_posterior(opt, img, burninfrac=0.5, rownum=rownum, colnum=colnum, nbins=100)
ax = gcf().axes
gcf().text(0.02, 0.9, "b.", fontsize=14, color="red")
ax[3].plot(img.x, img.f[rownum,:], "--y", alpha=0.8)
ax[5].plot(img.f[:,colnum], img.y,"--y", alpha=0.8)
savefig("post_s.png", dpi=300)
## plot comparisons
f, ax = plt.subplots(2, 2, sharex=true, sharey=true, figsize=(6.91, 6.94))
ax[1].imshow(img.f, extent=[img.x[1],img.x[end],img.y[end],img.y[1]])
ax[1].text(100, 170, "a. True image", color="w", fontsize=10, alpha=0.8)
ax[3].imshow(img.f, extent=[img.x[1],img.x[end],img.y[end],img.y[1]], alpha=0.0)
ax[3].scatter(Xtrain[1,:], Xtrain[2,:], c=ftrain, s=2)
ax[3].text(100, 170, "b. Noisy data", color="k", fontsize=10, alpha=0.8)
ax[3].text(100, 170, "b. Noisy data", color="w", fontsize=10, alpha=0.5)
ax[2].imshow(m, extent=[img.x[1],img.x[end],img.y[end],img.y[1]])
ax[2].text(100, 170, "c. Fixed λ", color="w", fontsize=10, alpha=0.8)
ax[4].imshow(mns, extent=[img.x[1],img.x[end],img.y[end],img.y[1]])
ax[4].text(100, 170, "d. Variable λ", color="w", fontsize=10, alpha=0.8)
for a in ax
   a.axis("off")
end
# f.tight_layout()
f.subplots_adjust(wspace=0, hspace=0)
savefig("compare_ns_s.png", dpi=300)
## plot stats
GeophysOperator.getchi2forall(opt, fsize=8, alpha=0.5)
gcf().text(0.02, 0.9, "a.", fontsize=14, color="red")
linidx = .!isnan.(img.d)
r = vec(img.d[linidx] - img.f[linidx])
χ² = length(r)#r'*r/img.σ^2
figure(1)
ax = gcf().axes
ax[2].set_ylim(χ²/2 - 100, χ²/2 + 300)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
figure(2)
ax = gcf().axes
ax[2].set_ylim(χ²/2 - 100, χ²/2 + 300)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
savefig("img_conv_s.png", dpi=300)
