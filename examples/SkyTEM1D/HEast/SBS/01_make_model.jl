srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, DelimitedFiles, Random, Statistics, Revise,
      AEM_VMD_HMD, SkyTEM1DInversion
## model fixed parts, i.e., air
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 100
# Note that the receiver and transmitter need to be in layer 1
zstart = 0.0
extendfrac, dz = 1.06, 2.
nlayers = 40
fdataname = "Line_115651_dbdt_gates_rangeidx_40.mat"
ρbg = 10.
ntimesperdecade = 10
nfreqsperdecade = 5
aem, znall = SkyTEM1DInversion.makeoperator(fdataname,
                       zfixed = zfixed,
                       ρfixed = ρfixed,
                       zstart = zstart,
                       extendfrac = extendfrac,
                       dz = dz,
                       ρbg = ρbg,
                       nlayers = nlayers,
                       ntimesperdecade = ntimesperdecade,
                       nfreqsperdecade = nfreqsperdecade,
                       showgeomplot = false,
                       plotfield = false)
## make options
using GP, TransD_GP
nmin, nmax = 2, 40
K = GP.Mat32()
demean = true
fbounds = [-0.5 2.5]
sdpos = 0.05
sdprop = 0.05
λ, δ = [2], 0.1
opt, optdummy = SkyTEM1DInversion.make_tdgp_statmode_opt(znall = znall,
                    fileprefix = "sounding",
                    nmin = nmin,
                    nmax = nmax,
                    K = K,
                    demean = demean,
                    sdpos = sdpos,
                    sdprop = sdprop,
                    fbounds = fbounds,
                    λ = λ,
                    δ = δ,
                    )
