any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using ImageRegression, TransD_GP, MCMC_Driver, MPI, Distributed

img = ImageRegression.Img(
          filename         = "4.2.01.tiff",
          dx               = 10.0,
          fractrain        = 0.02,
          dec              = 2,
          gausskernelwidth = 7)

d, sd, ftrain, Xtrain =  ImageRegression.get_training_data(img,
                   sdmaxfrac = 0.05, 
                   ybreak = 1000,
                   takeevery = 4)

ImageRegression.plot_data(ftrain, Xtrain, img)

Xall = ImageRegression.get_all_prediction_points(img)

opt_in = TransD_GP.Options(
              nmin = 2,
              nmax = 20,
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
              fdataname         = "2Dtest_smooth"
              )

opt_EM_in  = MCMC_Driver.EMoptions(sd=sd)

ImageRegression.calc_simple_RMS(d, img, opt_in, opt_EM_in, sd)

# actual run of McMC
nsamples, nchains, Tmax = 4001, 4, 2.5
mgr = MPI.start_main_loop(MPI.MPI_TRANSPORT_ALL)

addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
@time MCMC_Driver.main(opt_in, d, Tmax, nsamples, opt_EM_in)

rmprocs(workers())
MPI.stop_main_loop(mgr)

ImageRegression.plot_last_target_model(img, opt_in)
