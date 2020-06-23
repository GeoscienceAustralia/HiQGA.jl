using Statistics, LinearAlgebra
burninfrac = 0.9
temperaturenum = 8
x_ft = GeophysOperator.assembleTat1(optlog10λ, :x_ftrain, burninfrac=burninfrac, temperaturenum=temperaturenum)
n = GeophysOperator.assembleTat1(optlog10λ, :nodes, burninfrac=burninfrac, temperaturenum=temperaturenum)
nmodels = length(n)
i = nmodels
ftest, = GP.GPfit(optlog10λ.K, x_ft[i][2,1:n[i]]', x_ft[i][1,1:n[i]]',
    optlog10λ.xall, optlog10λ.λ², optlog10λ.δ, nogetvars=true, demean=optlog10λ.demean, p=2,
    my=mean(opt.fbounds, dims=2))
GeophysOperator.plot_posterior(csem, optlog10λ, burninfrac=burninfrac,
            pdfnormalize=true, cmappdf="inferno_r", temperaturenum=temperaturenum)
##
ax[2].plot(ftest,opt.xall[:])
