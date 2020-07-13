srcdir = dirname(dirname(dirname(pwd())))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, TransD_GP, Distributed, MCMC_Driver, GP, Random, Revise, GeophysOperator
##
img =     Img(
          filename         = "4.2.01.png",
          dx               = 10.0,
          fractrain        = 0.02,
          dec              = 2,
          gausskernelwidth = 7)
##
img.σ, ftrain, Xtrain =  get_training_data(img,
                   sdmaxfrac = 0.05,
                   ybreak = 1000,
                   takeevery = 4)

plot_data(ftrain, Xtrain, img)
calc_simple_RMS(img)
