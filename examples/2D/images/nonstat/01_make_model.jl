using PyPlot, Distributed, Random
using HiQGA.transD_GP
##
img =     transD_GP.Img(
          filename         = "4.2.01.png",
          dx               = 10.0,
          fractrain        = 0.02,
          dec              = 2,
          gausskernelwidth = 7)
##
transD_GP.get_image_data(img,
                   sdmaxfrac = 0.05,
                   ybreak = 1000,
                   takeevery = 4)

transD_GP.plot_image_data(img)
transD_GP.calc_image_RMS(img)
