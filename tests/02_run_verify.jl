srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Revise, CSEM1DEr, DelimitedFiles
##
rMin = 100 #m
rMax = 20000  #m
nRx = 150
freqs   = [0.1, 0.3, 0.7] #Hz
#receiver geometry
rRx   = collect(LinRange(rMin,rMax,nRx))  # Ranges to receivers   (m)
zRx = 1000.     # Depth of Receivers, if this is an array it calls normal FHT without lagged convolution
zTx = 975.
RxAzim=40.
TxDip=20.
# Note that the receiver depth needs to be in same model layer as transmitter.
# model
zfixed   = [-1e6,    0,      1000,   2000, 2100]
rhofixed = [1e13,    0.3,    1,      100,  1]
##
F = CSEM1DEr.RadialEr(zTx=[zTx],
                      rRx=rRx,
                      freqs=freqs,
                      zRx=[zRx])

##
figure()
step(log10.(rhofixed[2:end]), zfixed[2:end])
xlim(-1, 2)
ylim(1000, 3700)
grid()
##
krho = 2.0593e-10
rxno = 1
txno = 1
Ek = CSEM1DEr.getCSEM1DKernelsEr(F, krho, freqs[1], zfixed, rhofixed, rxno, txno)
##
@time Er = CSEM1Dkernels.getCSEM1DanisoHED(freqs, rRx, zRx, zTx, zfixed, [rhofixed rhofixed], 0,
                                            RxAzim = RxAzim, TxDip=TxDip)
figure()
semilogy(rRx, abs.(Er), label="julia")
Erm = readdlm("can_test_mat.txt",',',Complex{Float64})
semilogy(rRx, abs.(Erm), "--", label="mat")
