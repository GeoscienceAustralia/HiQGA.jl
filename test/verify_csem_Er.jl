using PyPlot, HiQGA.transD_GP, DelimitedFiles, Random, BenchmarkTools
## set geometry
rMin = 100 #m
rMax = 20000  #m
nRx = 150
freqs   = [0.1, 0.3, 0.7] #Hz
# receiver geometry
rRx   = collect(LinRange(rMin,rMax,nRx))  # Ranges to receivers   (m)
zRx = 975.     # Depth of Receivers, if this is an array it calls normal FHT without lagged convolution
zTx = 1000.
RxAzim=40.
TxDip=20.
# Note that the receiver depth needs to be in same model layer as transmitter.
# model
zfixed   = [-1e6,    0,      1000,   2000, 2100]
rhofixed = [1e13,    0.3,    1,      100,  1]
## set operator with lagged convolution
F = transD_GP.CSEM1DEr.RadialErLagged(zTx    = [zTx],
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = [zRx],
                      RxAzim = RxAzim,
                      TxDip  = TxDip)

## model one wavenumber
krho = 2.0593e-10
rxno = 1
txno = 1
Ek = transD_GP.CSEM1DEr.getCSEM1DKernelsEr!(F, krho, freqs[1], zfixed, rhofixed, rxno, txno)
## timing
@info "Lagged convolution time"
transD_GP.CSEM1DEr.getfield!(F, zfixed, rhofixed)
@btime transD_GP.CSEM1DEr.getfield!($F, $zfixed, $rhofixed)
# Discrete Hankel Tx
G = transD_GP.CSEM1DEr.RadialErDHT(rRx  = F.rRx,
                          freqs  = F.freqs,
                          zRx    = F.zRx[1]*ones(length(F.rRx)),
                          RxAzim = RxAzim*ones(length(F.rRx)),
                          zTx    = F.zTx[1]*ones(length(F.rRx)),
                          TxDip  = TxDip*ones(length(F.rRx))
                          )
@info "Vanilla DHT time"
transD_GP.CSEM1DEr.getfield!(G, zfixed, rhofixed)
@btime transD_GP.CSEM1DEr.getfield!($G, $zfixed, $rhofixed)
## plot
figure()
subplot(121)
step(log10.(rhofixed[2:end]), zfixed[2:end])
xlim(-1, 2)
ylim(1000, 3700)
grid()
gca().invert_yaxis()
subplot(122)
semilogy(rRx, abs.(F.Er), label="julia lagged HT")
semilogy(rRx, abs.(G.Er), label="julia vanilla DHT")
include("csem_Er_response.jl")
semilogy(rRx, abs.(Erm), "--", label="matlab 2012 code")
grid()
legend()
## random, many layers
Random.seed!(435)
nlayers = 100
z = [-1e6, 0, 1000, 1000 .+ cumsum(50*rand(nlayers-3))...]
ρ = [1e13, 0.3, 1, 10 .^(1 .+ 1.5*rand(nlayers-3))...]
figure()
subplot(121)
step(log10.(ρ[2:end]), z[2:end])
xlim(1,2.5 )
grid(); gca().invert_yaxis()
@info "timing for $nlayers layers"
transD_GP.CSEM1DEr.getfield!(F, z, ρ)
@btime transD_GP.CSEM1DEr.getfield!($F, $z, $ρ)
subplot(122)
semilogy(rRx, abs.(F.Er), label="julia")
