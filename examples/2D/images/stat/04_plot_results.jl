# convergence stats
close("all")
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
ax = gcf().axes;
r = vec(img.d[img.select] - img.f[img.select])
χ² = length(r)
figure(1)
ax[2].set_ylim(χ²/2 - 100, χ²/2 + 300)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
figure(2)
ax = gcf().axes
ax[2].set_ylim(χ²/2 - 100, χ²/2 + 300)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
## posterior marginals
rownum, colnum = 195, 85
transD_GP.plot_image_posterior(opt, img, burninfrac=0.5, rownum=rownum, colnum=colnum, nbins=100)
savefig("post_s.png", dpi=300)
