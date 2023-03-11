using PyPlot, Distributed, Random
using Revise, HiQGA.transD_GP
##
img =     transD_GP.Img(
          filename         = "../stationary/4.2.01.png",
          dx               = 10.0, # pixel spacing
          fractrain        = 0.02, # fraction to use as training data
          dec              = 2, # decimate original image by
          gausskernelwidth = 7 # blur to apply to image
          )
##
transD_GP.get_image_data(img,
                   sdmaxfrac = 0.05, # noise percentage
                   ybreak = 1000, # make sudden break at
                   takeevery = 4, # another decimnation to dataset after break
                   )

transD_GP.plot_image_data(img)
transD_GP.calc_image_RMS(img)
