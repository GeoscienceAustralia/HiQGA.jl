using PyPlot, Revise, transD_GP, DelimitedFiles, Random
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
## set operator
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
ntimes = 1000
t = time()
for i = 1:ntimes
    transD_GP.CSEM1DEr.getfield!(F, zfixed, rhofixed)
end
t = time() - t
@info "timing is $(t/ntimes) s"
## plot
figure()
subplot(121)
step(log10.(rhofixed[2:end]), zfixed[2:end])
xlim(-1, 2)
ylim(1000, 3700)
grid()
gca().invert_yaxis()
subplot(122)
semilogy(rRx, abs.(F.Er), label="julia")
# Erm = readdlm("can_test_mat.txt",',',Complex{Float64})
include("csem_Er_test.jl")
semilogy(rRx, abs.(Erm), "--", label="mat")
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
t = time()
for i = 1:ntimes
    transD_GP.CSEM1DEr.getfield!(F, z, ρ)
end
t = time() - t
@info "timing is $(t/ntimes) s"
subplot(122)
semilogy(rRx, abs.(F.Er), label="julia")
## Discrete Hankel Tx
G = transD_GP.CSEM1DEr.RadialErDHT(rRx  = F.rRx,
                          freqs  = F.freqs,
                          zRx    = F.zRx[1]*ones(length(F.rRx)),
                          RxAzim = RxAzim*ones(length(F.rRx)),
                          zTx    = F.zTx[1]*ones(length(F.rRx)),
                          TxDip  = TxDip*ones(length(F.rRx))
                          )
