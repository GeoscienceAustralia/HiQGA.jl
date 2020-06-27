srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Revise, AEM_VMD_HMD, Random, Statistics
## model
zfixed   = [-1e5,   0,    20,  50]
rho      = [1e12,   10,   1,   100]
nmax = 100
##  geometry
rRx = 17.
zRx = -37.0
zTx = -35.0
# Note that the receiver depth needs to be in same model layer as transmitter.
## times and ramp
HM_times = [
7.539e-05    9.639e-05   0.00012239   0.00015439   0.00019639   0.00024739   0.00031239   0.00039439   0.00049739   0.00062739   0.00079039   0.00099639    0.0012554    0.0015814    0.0019914    0.0025084    0.0031584    0.0039774    0.0050084    0.0063064    0.0079394
9.6e-05     0.000122     0.000154     0.000196     0.000247     0.000312     0.000394     0.000497     0.000627      0.00079     0.000996     0.001255     0.001581     0.001991     0.002508     0.003158     0.003977     0.005008     0.006306     0.007939     0.009739]'
HM_times = exp.(mean(log.(HM_times),dims=2))[:]
HM_ramp = [
-0.01    -0.008386     -0.00638    -0.003783            0     3.96e-07    7.782e-07    1.212e-06     3.44e-06    1.981e-05    3.619e-05    3.664e-05    3.719e-05    3.798e-05    3.997e-05         0.01
0       0.4568       0.7526       0.9204            1       0.9984       0.9914       0.9799       0.9175       0.4587     0.007675     0.003072    0.0008319     0.000119            0            0]'
lowpassfcs = [300000, 450000.0]
ntimesperdecade = 10
nfreqsperdecade = 10
##
F = AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      lowpassfcs = lowpassfcs,
                      times  = HM_times,
                      ramp   = HM_ramp,
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      zRx    = zRx)

##
@time AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
##
figure()
loglog(F.times,4*pi*1e-7*abs.(F.dBzdt))
# xlim(extrema(F.freqs))
# legend()
grid()
## timing
ntimes = 1000
t = time()
for i = 1:ntimes
    AEM_VMD_HMD.getfieldFD!(F, zfixed, rho)
end
t = time() - t
@info "timing is $(t/ntimes) s"
