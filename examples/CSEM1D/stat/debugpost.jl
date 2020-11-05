using Statistics, LinearAlgebra
burninfrac = 0.9
temperaturenum = 6
x_ft = transD_GP.assembleTat1(opt, :x_ftrain, burninfrac=burninfrac, temperaturenum=temperaturenum)
n = transD_GP.assembleTat1(opt, :nodes, burninfrac=burninfrac, temperaturenum=temperaturenum)
nmodels = length(n)
i = nmodels
ftest, = transD_GP.GP.GPfit(opt.K, x_ft[i][2,1:n[i]]', x_ft[i][1,1:n[i]]',
    opt.xall, opt.λ², opt.δ, nogetvars=true, demean=opt.demean, p=2,
    my=mean(opt.fbounds, dims=2))
transD_GP.plot_posterior(csem, opt, burninfrac=burninfrac,
            pdfnormalize=true, cmappdf="inferno_r", temperaturenum=temperaturenum)
##
ax[2].plot(ftest,opt.xall[:])
