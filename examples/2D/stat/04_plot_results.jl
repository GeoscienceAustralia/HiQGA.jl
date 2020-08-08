GeophysOperator.plot_image_data(ftrain, Xtrain, img)
savefig("image2d.png", dpi=300)
GeophysOperator.plot_image_posterior(opt, img, burninfrac=0.5, rownum=168, colnum=60, nbins=100)
savefig("post_s.png", dpi=300)
