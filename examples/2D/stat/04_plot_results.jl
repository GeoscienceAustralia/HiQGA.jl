GeophysOperator.plot_image_data(ftrain, Xtrain, img)
savefig("image2d.png", dpi=300)
rownum, colnum = 168, 60
GeophysOperator.plot_image_posterior(opt, img, burninfrac=0.5, rownum=rownum, colnum=colnum, nbins=100)
savefig("post_s.png", dpi=300)
