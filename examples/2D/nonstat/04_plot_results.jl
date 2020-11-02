## get nonstat stuff first
close("all")
using Statistics
M = transD_GP.assembleTat1(opt, :fstar, temperaturenum=1)
mns = reshape(mean(M), length(img.y), length(img.x))
stdns = reshape(std(M), length(img.y), length(img.x))
##
rownum, colnum = 195, 85
transD_GP.plot_image_posterior(opt, optlog10λ, img, burninfrac=0.5, rownum=rownum, colnum=colnum, nbins=100)
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
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
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
transD_GP.getchi2forall(optlog10λ, fsize=8, alpha=0.5)
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
