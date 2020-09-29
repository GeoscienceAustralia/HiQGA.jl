srcdir = dirname(dirname(dirname(pwd())))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, DelimitedFiles, Random, Statistics, Revise,
      AEM_VMD_HMD, SkyTEM1DInversion, GeophysOperator
## model fixed parts, i.e., air
Random.seed!(23)
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 100
# Note that the receiver and transmitter need to be in layer 1
zstart = 0.0
extendfrac, dz = 1.06, 2.
zall, znall, zboundaries = GeophysOperator.setupz(zstart, extendfrac, dz=dz, n=40)
z, ρ, nfixed = GeophysOperator.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
ρ[z.>=zstart] .= 10.
##  geometry and modeling parameters
rRx = 13.05
zRx = -40.78
zTx = -37.9
rTx = 10.4184
lowpassfcs = [300000, 155000]
ntimesperdecade = 10
nfreqsperdecade = 5
# Note that the receiver depth needs to be in same model layer as transmitter.
## LM times and ramp
LM_times = [0.00001463	0.0000182;
			0.00001863	0.0000232;
			0.00002363	0.0000292;
			0.00002963	0.0000372;
			0.00003763	0.0000472;
			0.00004763	0.0000602;
			0.00006063	0.0000762;
			0.00007663	0.0000962;
			0.00009663	0.0001212;
			0.00012163	0.0001522;
			0.00015263	0.0001922;
			0.00019263	0.0002432;
			0.00024363	0.0003062;
			0.00030663	0.0003872;
			0.00038763	0.0004882;
			0.00048863	0.0006152;
			0.00061563	0.0007762;
			0.00077663	0.0009782]
LM_times = exp.(mean(log.(LM_times),dims=2))[:]
LM_ramp =[-8.0000E-04	0.0000E+00;
			-7.1847E-04	1.2266E-01;
			-6.5990E-04	1.9609E-01;
			-3.8062E-04	5.2228E-01;
			-1.6838E-04	7.7866E-01;
			-5.1442E-05	9.3110E-01;
			0.0000E+00	1.0000E+00;
			1.2800E-07	9.9627E-01;
			6.2400E-07	9.7994E-01;
			1.1840E-06	9.3093E-01;
			2.9600E-06	6.7610E-01;
			5.1200E-06	3.5920E-01;
			6.2400E-06	2.2688E-01;
			7.5200E-06	1.2234E-01;
			8.8000E-06	5.6996E-02;
			1.0048E-05	2.2692E-02;
			1.2000E-05	4.7230E-03;
			1.2800E-05	0.0000E+00;
			1.0182E-03	0.0000E+00]
## HM times and ramp
HM_times = [4.26630E-04	4.46200E-04;
			4.46630E-04	4.71200E-04;
			4.71630E-04	5.02200E-04;
			5.02630E-04	5.42200E-04;
			5.42630E-04	5.93200E-04;
			5.93630E-04	6.56200E-04;
			6.56630E-04	7.37200E-04;
			7.37630E-04	8.38200E-04;
			8.38630E-04	9.65200E-04;
			9.65630E-04	1.12620E-03;
			1.12663E-03	1.32820E-03;
			1.32863E-03	1.58320E-03;
			1.58363E-03	1.90520E-03;
			1.90563E-03	2.31120E-03;
			2.31163E-03	2.82220E-03;
			2.82263E-03	3.46820E-03;
			3.46863E-03	4.26020E-03;
			4.26063E-03	5.22820E-03;
			5.22863E-03	6.41320E-03;
			6.41363E-03	7.86520E-03;
			7.86563E-03	9.64120E-03;
			9.64163E-03	1.18212E-02;
			1.18216E-02	1.44912E-02]
HM_times = exp.(mean(log.(HM_times),dims=2))[:]
HM_ramp = [-5.00000E-03 0.00000E+00;
    -4.79532E-03 6.52344E-01;
    -4.69386E-03 9.33594E-01;
	-3.82897E-03 9.49219E-01;
	-2.13784E-03 9.77344E-01;
	-7.26966E-04 9.94531E-01;
	-2.19939E-06 1.00000E+00;
	1.04316E-05 9.81961E-01;
	5.79841E-05 8.39282E-01;
	1.48063E-04 5.57539E-01;
	2.31183E-04 2.88357E-01;
	2.79509E-04 1.29935E-01;
	3.15850E-04 1.35143E-02;
	3.28993E-04 0.00000E+00;
	1.50000E-02 0.00000E+00]
## LM operator
Flm = AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      lowpassfcs = lowpassfcs,
                      times  = LM_times,
                      ramp   = LM_ramp,
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx)
## HM operator
Fhm = AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
                      lowpassfcs = lowpassfcs,
                      times  = HM_times,
                      ramp   = HM_ramp,
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx)
## data and high altitude noise
LM_data = 1e-12*[7.202820E+03   4.943799E+03   3.077346E+03   1.932652E+03   1.209026E+03   7.387355E+02   4.463695E+02   2.659717E+02   1.571060E+02   9.194739E+01   5.221553E+01   2.882322E+01   1.576041E+01   8.561270E+00   4.659970E+00   2.554020E+00   1.399650E+00   7.654700E-01]'
HM_data = 1e-12*[3.324394E+01   2.522483E+01   1.890557E+01   1.373576E+01   9.740710E+00   6.743540E+00   4.606930E+00   3.036280E+00   1.992620E+00   1.262880E+00   7.893800E-01   4.761200E-01   2.914300E-01   1.748700E-01   9.983000E-02   6.034000E-02   3.520000E-02   2.008000E-02   1.102000E-02   4.940000E-03   1.510000E-03  -1.300000E-04  -3.900000E-04]'
LM_noise = 1e-12*[7.396625e+00    4.862768e+00    3.903425e+00    2.763621e+00    2.266908e+00    1.737678e+00    1.361518e+00    1.081970e+00    8.484297e-01    6.708301e-01    5.049530e-01    3.889236e-01    3.091977e-01    2.402437e-01    1.942889e-01    1.556078e-01    1.346912e-01    1.260919e-01]'
HM_noise = 1e-12*[6.044069e-02    3.294512e-02    2.174927e-02    1.393709e-02    1.029002e-02    8.436072e-03    6.701908e-03    5.082382e-03    4.312096e-03    3.399737e-03    2.691244e-03    2.320745e-03    2.018198e-03    1.759839e-03    1.586694e-03    1.435991e-03    1.282377e-03    1.097777e-03    9.950349e-04    8.986489e-04    8.252542e-04    7.411206e-04    6.240619e-04]'
## create operator
dlow, dhigh, σlow, σhigh = (LM_data, HM_data, 0.03*LM_data+LM_noise, 0.03*HM_data+HM_noise)./SkyTEM1DInversion.μ₀
aem = GeophysOperator.dBzdt(Flm, Fhm, vec(dlow), vec(dhigh),
                                  vec(σlow), vec(σhigh), z=z, ρ=ρ, nfixed=nfixed)
GeophysOperator.SkyTEM1DInversion.plotmodelfield!(Flm, Fhm, z, ρ, dlow, dhigh, σlow, σhigh;
                      figsize=(12,4), nfixed=nfixed, dz=dz, extendfrac=extendfrac)
