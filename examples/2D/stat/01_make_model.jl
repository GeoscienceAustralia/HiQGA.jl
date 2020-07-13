srcdir = dirname(dirname(dirname(pwd())))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Distributed, Random
using Revise, GeophysOperator, TransD_GP, MCMC_Driver, GP
##
img =     GeophysOperator.Img(
          filename         = "4.2.01.png",
          dx               = 10.0,
          fractrain        = 0.02,
          dec              = 2,
          gausskernelwidth = 7)
##
img.σ, ftrain, Xtrain = GeophysOperator.get_image_data(img,
                   sdmaxfrac = 0.05,
                   ybreak = 1000,
                   takeevery = 4)

GeophysOperator.plot_image_data(ftrain, Xtrain, img)
GeophysOperator.calc_image_RMS(img)
