any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using TransD_GP, Distributed, Revise, GeophysOperator, MCMC_Driver
##
img =     Img(
          filename         = "4.2.01.tiff",
          dx               = 10.0,
          fractrain        = 0.02,
          dec              = 2,
          gausskernelwidth = 7)
##
sd, ftrain, Xtrain =  get_training_data(img,
                   sdmaxfrac = 0.05,
                   ybreak = 1000,
                   takeevery = 4)

plot_data(ftrain, Xtrain, img)

Xall = get_all_prediction_points(img)

opt_in = TransD_GP.Options(
              nmin = 2,
              nmax = 75,
              xbounds           = [img.x[1] img.x[end];img.y[1] img.y[end]],
              fbounds           = [-1 2],
              xall              = Xall,
              λ                 = [150.0, 150.0],
              δ                 = 0.2,
              demean            = true,
              save_freq         = 500,
              dispstatstoscreen = false,
              sdev_prop         = 0.1,
              sdev_pos          = [10.0, 10.0],
              pnorm             = 2,
              debug             = false,
              fdataname         = "2Dtest_smooth",
              quasimultid       = false
              )

calc_simple_RMS(img, sd)

# actual run of McMC
nsamples, nchains, nchainsatone = 501, 8, 1
Tmax = 2.5

addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
##
@time MCMC_Driver.main(opt_in, img, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
##
rmprocs(workers())
## plot
#ImageRegression.getchi2forall(opt_in)
#ImageRegression.plot_last_target_model(img, opt_in)
