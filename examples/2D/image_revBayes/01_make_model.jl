using PyPlot, Distributed, Random
using HiQGA.transD_GP
##
img =     transD_GP.Img(
          filename         = "bayesSmall.png",
          dx               = 10.0, # define spacing
          dec              = 5, # decimate image by
          fractrain        = 0.005, # fraction of pixels to train
          gausskernelwidth = 6)
##
transD_GP.get_image_data(img,
                   sdmaxfrac = 0.05, # noise fraction
                   ybreak = 400, # make a break in the image at
                   rseed = 1, # random seed
                   takeevery = 4 # decimate after the break by
                   )

transD_GP.plot_image_data(img)
transD_GP.calc_image_RMS(img)
