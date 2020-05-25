module CSEM1Dkernels
  using Dierckx, LinearAlgebra
  export stacks, getCurlyR, getCSEM1DKernelsAnisoHED, getCSEM1DanisoHED

  function stacks(pz,iTxLayer,r,z,omega)

    #The last and first layer thicknesses are infinite

    dz          = diff(z)
    nz          = length(z)
    d           = zeros(1,nz)
    d[1]        = 1e60
    d[2:nz-1]   = dz[2:end]
    d[nz]       = 1e60

  # Capital R is for a stack

  # Starting from the bottom up, for Rs_down
    Rlowerstack=0. *im;
    for k = (nz-1):-1:iTxLayer
        Rs_d = (r[k] + Rlowerstack*exp(2im*omega*pz[k+1]*d[k+1])) /
              (1. + r[k]*Rlowerstack *
              exp(2im*omega*pz[k+1]*d[k+1]))
        Rlowerstack=Rs_d
    end

    Rupperstack=0. *im;
  # Starting from the top down for Rs_up
    for k   = 2:iTxLayer
        Rs_u = (-r[k-1] + Rupperstack*exp(2im*omega*pz[k-1]*d[k-1])) /
              (1. - r[k-1]*Rupperstack *
              exp(2im*omega*pz[k-1]*d[k-1]))
        Rupperstack=Rs_u
    end

    return Rupperstack,Rlowerstack

  end

  function getCurlyR(Rs_u,Rs_d,pz,zR,z,iTxLayer,omega)
    d=z[iTxLayer+1]-z[iTxLayer]
    if (zR>=0)
        finRA = (1. + exp(-im*omega*pz*2*z[iTxLayer])*Rs_u) *
                (exp(im*omega*pz*zR) + exp( im*omega*pz*(2*z[iTxLayer+1] - zR))*Rs_d) /
                (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))

        finRD = (1. - exp(-im*omega*pz*2*z[iTxLayer])*(-Rs_u)) *
                    (exp(im*omega*pz*zR) + exp( im*omega*pz*(2*z[iTxLayer+1] - zR))*(-Rs_d)) /
                    (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))

        finRB = (1. - exp(-im*omega*pz*2*z[iTxLayer])*Rs_u) *
                (exp(im*omega*pz*zR) + exp( im*omega*pz*(2*z[iTxLayer+1] - zR))*Rs_d) /
                (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))
    else
        finRA = (1. + exp(im*omega*pz*2*z[iTxLayer+1])*Rs_d) *
                (exp(-im*omega*pz*zR) + exp(-im*omega*pz*(2*z[iTxLayer] - zR))*Rs_u) /
                (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))

        finRD = (-1. + exp(im*omega*pz*2*z[iTxLayer+1])*(-Rs_d)) *
                (exp(-im*omega*pz*zR) + exp(-im*omega*pz*(2*z[iTxLayer] - zR))*(-Rs_u)) /
                (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))

        finRB = (-1. + exp(1im*omega*pz*2*z[iTxLayer+1])*Rs_d) *
                (exp(-im*omega*pz*zR) + exp(-im*omega*pz*(2*z[iTxLayer] - zR))*Rs_u) /
                (1. - Rs_u*Rs_d*exp(2im*omega*pz*(d)))
    end

    return finRA, finRD, finRB
  end

  function getCSEM1DKernelsAnisoHED(krho,f,zRx,zTx,z,rho)

    z       = z.-zTx
    zRx     = zRx-zTx
    zTx     = 0.;
    # will need these variables
    mu      = 4 *pi*1e-7
    eps0    = 8.854e-12
    omega   = 2. *pi*f
    epsc    = eps0 .+ 1im./(rho*omega)

    pzI  = 1im*zeros(size(epsc,1))
    pzII = similar(pzI)
    zRxFlag = Array{Int}(undef, length(z))
    zTxFlag = similar(zRxFlag)
    for layer in eachindex(z)
      pzI[layer]  = sqrt(mu*epsc[layer,1] - (krho/omega)^2)
      pzII[layer]        = sqrt(mu*epsc[layer,1] - epsc[layer,1]/epsc[layer,2]*(krho/omega)^2)
      # wavenubmer sanity check
      imag(pzI[layer])  < 0.0 && ( pzI[layer]*=-1.)
      imag(pzII[layer]) < 0.0 && ( pzII[layer]*=-1.)
      #where are the Tx and Rx
      zRxFlag[layer] =  round(Int,z[layer]<zTx)
      zTxFlag[layer] =  round(Int,z[layer]<zRx)
    end
    # reflection coefficients (downward) for an intfc: pzI for TE, pzII for TM
    rTE = 1im*zeros(length(pzI)-1)
    rTM = (similar(rTE))
    for intfc in 1:length(pzI)-1
      rTE[intfc] = (pzI[intfc] - pzI[intfc+1])/(pzI[intfc] + pzI[intfc+1])
      rTM[intfc] = (epsc[intfc,1]*pzII[intfc+1] - epsc[intfc+1,1]*pzII[intfc]) /
                   (epsc[intfc,1]*pzII[intfc+1] + epsc[intfc+1,1]*pzII[intfc])
    end
    # Find the layer containing the transmitter:
    iTxLayer = sum(zTxFlag)
    if iTxLayer ==  0
      iTxLayer = length(z) # in the last layer
    end

    # Find the layer containing the receiver:
    iRxLayer = sum(zRxFlag)
    if iRxLayer == 0
      iRxLayer = length(z) # in the last layer
    end

    @assert iTxLayer == iRxLayer "Error, transmitter and receiver need to be located in the same layer! Stopping"

    # TE mode, pzI
    Rs_u,Rs_d       = stacks(pzI,iTxLayer,rTE,z,omega)
    curlyRA,curlyRD = getCurlyR(Rs_u,Rs_d,pzI[iTxLayer],zRx,z,iTxLayer,omega)
    gA_TE             = mu/pzI[iTxLayer]*curlyRA
    gD_TE             = curlyRD

    # TM mode, pzII
    Rs_u,Rs_d       = stacks(pzII,iTxLayer,rTM,z,omega)
    curlyRA,curlyRD,curlyRB = getCurlyR(Rs_u,Rs_d,pzII[iTxLayer],zRx,z,iTxLayer,omega)
    gA_TM             = pzII[iTxLayer]/epsc[iTxLayer,1]*curlyRA
    gB_TM             = curlyRB
    gD_TM             = curlyRD

    # Kernels according to Loseth, without the bessel functions multiplied
    #Er from HED and VED
    J0                = -krho*gA_TM/4/pi
    J1                = -(gA_TE - gA_TM)/4/pi
    J1V               = krho^2*gB_TM*1im/(4*pi*omega*epsc[iTxLayer,2]);
    ErKernels         = [J0; J1; J1V]

     #Eb
    J0                = krho*gA_TE/4/pi
    J1                = -(gA_TE - gA_TM)/4/pi
    EbKernels         = [J0; J1]

    #Hr
    J0                = -krho*gD_TE/4/pi
    J1                = (gD_TE - gD_TM)/4/pi
    HrKernels         = [J0; J1]

    #Hb
    J0                = -krho*gD_TM/4/pi
    J1                = -(gD_TE - gD_TM)/4/pi
    HbKernels         = [J0; J1]

    #Ez
    EzKernels         = krho^2*gD_TM*im/(4*pi*omega*epsc[iTxLayer,2]);
    # needs multiplication by 1i/omega/epsc[2,iTxLayer] after hankel
    # transform (done here)

    #Hz
    HzKernels         = krho^2*gA_TE*im/(omega*mu*4*pi);
    # needs multiplication by 1i/omega/mu after hankel
    # transform (done here)

    return ErKernels, EbKernels, HrKernels, HbKernels, EzKernels, HzKernels
  end

  #Digital filter weights for Hankel transform
  const Filter_J0 = [   1.5020099e-03  -1.0381698e-02   3.6840860e-02  -8.9903380e-02   1.7082287e-01  -2.7115750e-01   3.7649328e-01  -4.7220779e-01   5.4778211e-01  -5.9823517e-01   6.2345792e-01  -6.2650437e-01   6.1197225e-01  -5.8470174e-01   5.4911056e-01  -5.0878862e-01   4.6652106e-01  -4.2426479e-01   3.8339539e-01  -3.4472345e-01   3.0876891e-01  -2.7570439e-01   2.4562331e-01  -2.1839207e-01   1.9393194e-01  -1.7198163e-01   1.5242270e-01  -1.3494946e-01   1.1946763e-01  -1.0565871e-01   9.3482356e-02  -8.2612525e-02   7.3077764e-02  -6.4536481e-02   5.7096311e-02  -5.0384860e-02   4.4599557e-02  -3.9316995e-02   3.4838545e-02  -3.0664946e-02   2.7221072e-02  -2.3901586e-02   2.1281567e-02  -1.8612485e-02   1.6655386e-02  -1.4472420e-02   1.3057748e-02  -1.1226286e-02   1.0266819e-02  -8.6738022e-03   8.1103369e-03  -6.6573628e-03   6.4551190e-03  -5.0524015e-03   5.1989014e-03  -3.7597430e-03   4.2640545e-03  -2.6994966e-03   3.5928001e-03  -1.8061291e-03   3.1436471e-03  -1.0244200e-03   2.8888298e-03  -3.0605071e-04   2.8125907e-03   3.9337862e-04   2.9102041e-03   1.1170885e-03   3.1876763e-03   1.9097806e-03   3.6621020e-03   2.8203784e-03   4.3626886e-03   3.9050025e-03   5.3324923e-03   5.2303399e-03   6.6309348e-03   6.8775541e-03   8.3371690e-03   8.9468565e-03   1.0554327e-02   1.1562766e-02   1.3414536e-02   1.4879837e-02   1.7084241e-02   1.9088093e-02   2.1768514e-02   2.4416121e-02   2.7711234e-02   3.1127164e-02   3.5184140e-02   3.9497913e-02   4.4449762e-02   4.9758433e-02   5.5667457e-02   6.1949790e-02   6.8682065e-02   7.5616744e-02   8.2584984e-02   8.9193635e-02   9.4874687e-02   9.8891682e-02   1.0004295e-01   9.7016844e-02   8.7789596e-02   7.0509768e-02   4.2778853e-02   3.5584533e-03  -4.7210453e-02  -1.0489788e-01  -1.6020950e-01  -1.9459782e-01  -1.8490775e-01  -1.0754165e-01   3.6037727e-02   1.9759013e-01   2.6132313e-01   1.1713997e-01  -1.8758779e-01  -3.0238115e-01   4.8163136e-02   3.6399530e-01  -1.4910233e-01  -2.6373490e-01   4.0362662e-01  -3.1409795e-01   1.8179369e-01  -9.0738718e-02   4.2946488e-02  -2.0586136e-02   1.0392667e-02  -5.6117848e-03   3.2402026e-03  -1.9858724e-03   1.2807317e-03  -8.6253792e-04   6.0296591e-04  -4.3548937e-04   3.2375892e-04  -2.4698212e-04   1.9279062e-04  -1.5357912e-04   1.2453788e-04  -1.0255126e-04   8.5558482e-05  -7.2170928e-05   6.1436863e-05  -5.2693349e-05   4.5471255e-05  -3.9433284e-05   3.4332971e-05  -2.9987212e-05   2.6257657e-05  -2.3037978e-05   2.0245071e-05  -1.7812926e-05   1.5688306e-05  -1.3827679e-05   1.2195005e-05  -1.0760111e-05   9.4974858e-06  -8.3853711e-06   7.4050613e-06  -6.5403683e-06   5.7772102e-06  -5.1032933e-06   4.5078641e-06  -3.9815112e-06   3.5159993e-06  -3.1041250e-06   2.7395854e-06  -2.4168588e-06   2.1310948e-06  -1.8780186e-06   1.6538488e-06  -1.4552319e-06   1.2791891e-06  -1.1230742e-06   9.8453631e-07  -8.6148574e-07   7.5206243e-07  -6.5460810e-07   5.6764452e-07  -4.8985885e-07   4.2009674e-07  -3.5736213e-07   3.0082216e-07  -2.4981510e-07   2.0385823e-07  -1.6265189e-07   1.2607417e-07  -9.4158418e-08   6.7043911e-08  -4.4891091e-08   2.7761326e-08  -1.5480404e-08   7.5327300e-09  -3.0524770e-09   9.5877856e-10  -2.0575286e-10   2.2414417e-11]'
  const Filter_J1 = [   4.7827871e-10  -2.9784176e-09   9.7723833e-09  -2.2382341e-08   4.0446774e-08  -6.1734816e-08   8.3293912e-08  -1.0249454e-07   1.1780780e-07  -1.2870061e-07   1.3559243e-07  -1.3921011e-07   1.4065746e-07  -1.4074882e-07   1.4051721e-07  -1.4040747e-07   1.4127886e-07  -1.4315596e-07   1.4689283e-07  -1.5210917e-07   1.5989802e-07  -1.6940918e-07   1.8227089e-07  -1.9695857e-07   2.1603952e-07  -2.3691321e-07   2.6369843e-07  -2.9202021e-07   3.2852445e-07  -3.6589095e-07   4.1487501e-07  -4.6327137e-07   5.2852697e-07  -5.9034984e-07   6.7710212e-07  -7.5512943e-07   8.7062473e-07  -9.6788531e-07   1.1222658e-06  -1.2417228e-06   1.4493467e-06  -1.5932457e-06   1.8747046e-06  -2.0433320e-06   2.4285695e-06  -2.6179926e-06   3.1511730e-06  -3.3492412e-06   4.0964187e-06  -4.2758276e-06   5.3371198e-06  -5.4435424e-06   6.9725767e-06  -6.9045563e-06   9.1397025e-06  -8.7148373e-06   1.2029591e-05  -1.0927977e-05   1.5912527e-05  -1.3582560e-05   2.1176227e-05  -1.6678206e-05   2.8384979e-05  -2.0132088e-05   3.8372045e-05  -2.3702184e-05   5.2385309e-05  -2.6854374e-05   7.2318582e-05  -2.8535362e-05   1.0108123e-04  -2.6788478e-05   1.4319185e-04  -1.8108424e-05   2.0573562e-04   3.6361649e-06   2.9991265e-04   4.8993332e-05   4.4354734e-04   1.3589102e-04   6.6515824e-04   2.9451992e-04   1.0105554e-03   5.7533965e-04   1.5535077e-03   1.0621193e-03   2.4128970e-03   1.8929698e-03   3.7800178e-03   3.2937343e-03   5.9612180e-03   5.6295935e-03   9.4422527e-03   9.4810228e-03   1.4979159e-02   1.5745093e-02   2.3708966e-02   2.5740591e-02   3.7232783e-02   4.1225009e-02   5.7507103e-02   6.4044643e-02   8.6091797e-02   9.4717140e-02   1.2172497e-01   1.2853597e-01   1.5450777e-01   1.4755964e-01   1.5621399e-01   1.1147621e-01   7.7489831e-02  -2.7628267e-02  -1.0198730e-01  -2.2039890e-01  -2.1185763e-01  -1.6052415e-01   9.1649026e-02   2.3792824e-01   2.6075778e-01  -1.5662188e-01  -2.8932082e-01   1.3148519e-02   4.2691303e-01  -4.0005050e-01   1.1513789e-01   9.3748244e-02  -1.6037231e-01   1.5071858e-01  -1.2120369e-01   9.4110656e-02  -7.3742238e-02   5.9038568e-02  -4.8288118e-02   4.0197054e-02  -3.3919788e-02   2.8918247e-02  -2.4845272e-02   2.1470450e-02  -1.8635828e-02   1.6229579e-02  -1.4170085e-02   1.2396084e-02  -1.0860414e-02   9.5259444e-03  -8.3628577e-03   7.3468030e-03  -6.4576043e-03   5.6783440e-03  -4.9946949e-03   4.3944258e-03  -3.8670264e-03   3.4034180e-03  -2.9957261e-03   2.6370977e-03  -2.3215540e-03   2.0438677e-03  -1.7994617e-03   1.5843227e-03  -1.3949289e-03   1.2281870e-03  -1.0813787e-03   9.5211407e-04  -8.3829103e-04   7.3806018e-04  -6.4979407e-04   5.7206023e-04  -5.0359765e-04   4.4329604e-04  -3.9017773e-04   3.4338167e-04  -3.0214942e-04   2.6581281e-04  -2.3378310e-04   2.0554144e-04  -1.8063040e-04   1.5864668e-04  -1.3923451e-04   1.2208003e-04  -1.0690622e-04   9.3468580e-05  -8.1551329e-05   7.0964175e-05  -6.1539592e-05   5.3130609e-05  -4.5609106e-05   3.8864649e-05  -3.2803856e-05   2.7350297e-05  -2.2444816e-05   1.8046076e-05  -1.4130827e-05   1.0693107e-05  -7.7412053e-06   5.2910576e-06  -3.3552268e-06   1.9282956e-06  -9.7253713e-07   4.1100808e-07  -1.3553176e-07   3.0748588e-08  -3.5668195e-09]'
  const Filter_base = [   4.1185887e-06   4.6623078e-06   5.2778064e-06   5.9745606e-06   6.7632974e-06   7.6561600e-06   8.6668947e-06   9.8110623e-06   1.1106278e-05   1.2572483e-05   1.4232251e-05   1.6111134e-05   1.8238059e-05   2.0645772e-05   2.3371342e-05   2.6456730e-05   2.9949438e-05   3.3903239e-05   3.8379005e-05   4.3445642e-05   4.9181157e-05   5.5673850e-05   6.3023682e-05   7.1343808e-05   8.0762323e-05   9.1424231e-05   1.0349368e-04   1.1715649e-04   1.3262301e-04   1.5013135e-04   1.6995107e-04   1.9238731e-04   2.1778548e-04   2.4653662e-04   2.7908337e-04   3.1592681e-04   3.5763416e-04   4.0484754e-04   4.5829384e-04   5.1879590e-04   5.8728520e-04   6.6481616e-04   7.5258245e-04   8.5193528e-04   9.6440425e-04   1.0917209e-03   1.2358454e-03   1.3989966e-03   1.5836864e-03   1.7927581e-03   2.0294306e-03   2.2973477e-03   2.6006340e-03   2.9439590e-03   3.3326083e-03   3.7725655e-03   4.2706040e-03   4.8343916e-03   5.4726080e-03   6.1950791e-03   7.0129278e-03   7.9387456e-03   8.9867860e-03   1.0173184e-02   1.1516206e-02   1.3036528e-02   1.4757557e-02   1.6705789e-02   1.8911218e-02   2.1407799e-02   2.4233968e-02   2.7433236e-02   3.1054859e-02   3.5154593e-02   3.9795557e-02   4.5049202e-02   5.0996412e-02   5.7728748e-02   6.5349859e-02   7.3977077e-02   8.3743226e-02   9.4798660e-02   1.0731359e-01   1.2148069e-01   1.3751806e-01   1.5567263e-01   1.7622389e-01   1.9948824e-01   2.2582385e-01   2.5563618e-01   2.8938422e-01   3.2758753e-01   3.7083428e-01   4.1979029e-01   4.7520927e-01   5.3794444e-01   6.0896164e-01   6.8935424e-01   7.8035994e-01   8.8337984e-01   1.0000000e+00   1.1320159e+00   1.2814599e+00   1.4506330e+00   1.6421396e+00   1.8589280e+00   2.1043360e+00   2.3821418e+00   2.6966223e+00   3.0526193e+00   3.4556135e+00   3.9118093e+00   4.4282302e+00   5.0128269e+00   5.6745996e+00   6.4237368e+00   7.2717720e+00   8.2317613e+00   9.3184844e+00   1.0548672e+01   1.1941264e+01   1.3517701e+01   1.5302252e+01   1.7322392e+01   1.9609223e+01   2.2197951e+01   2.5128433e+01   2.8445785e+01   3.2201080e+01   3.6452134e+01   4.1264394e+01   4.6711949e+01   5.2878668e+01   5.9859491e+01   6.7761894e+01   7.6707539e+01   8.6834152e+01   9.8297638e+01   1.1127449e+02   1.2596448e+02   1.4259380e+02   1.6141844e+02   1.8272824e+02   2.0685126e+02   2.3415891e+02   2.6507161e+02   3.0006526e+02   3.3967864e+02   3.8452161e+02   4.3528457e+02   4.9274904e+02   5.5779973e+02   6.3143815e+02   7.1479801e+02   8.0916269e+02   9.1598501e+02   1.0369096e+03   1.1737981e+03   1.3287581e+03   1.5041752e+03   1.7027502e+03   1.9275403e+03   2.1820062e+03   2.4700656e+03   2.7961535e+03   3.1652901e+03   3.5831587e+03   4.0561925e+03   4.5916743e+03   5.1978481e+03   5.8840466e+03   6.6608341e+03   7.5401699e+03   8.5355920e+03   9.6624257e+03   1.0938019e+04   1.2382011e+04   1.4016633e+04   1.5867051e+04   1.7961754e+04   2.0332991e+04   2.3017268e+04   2.6055913e+04   2.9495707e+04   3.3389608e+04   3.7797566e+04   4.2787445e+04   4.8436067e+04   5.4830397e+04   6.2068879e+04   7.0262956e+04   7.9538782e+04   9.0039163e+04   1.0192576e+05   1.1538158e+05   1.3061378e+05   1.4785687e+05   1.6737633e+05   1.8947266e+05   2.1448605e+05   2.4280162e+05]

function getCSEM1DanisoHED(freqs::Array{T,1}, rRx::Array{T,1}, zRx::T, zTx::T, z::Array{T,1}, rho::Array{T,2}) where T<:AbstractFloat
    Er = 1.0im*zeros(length(rRx),length(freqs)) #space domain
    Eb = 1.0im*zeros(length(rRx),length(freqs))
    Hr = 1.0im*zeros(length(rRx),length(freqs))
    Hb = 1.0im*zeros(length(rRx),length(freqs))
    Ez = 1.0im*zeros(length(rRx),length(freqs))
    Hz = 1.0im*zeros(length(rRx),length(freqs))

    #lagged convolution, good for all Rx depths same
    lambdaR = 1.
    rR = 1.
    rMax   = maximum(rRx)
    rMin   = minimum(rRx)
    lamMin = Filter_base[1]/rMax
    lamMax = Filter_base[end]/rMin
    filterSpacing = Filter_base[2] / Filter_base[1] # the filters are spaced by a multiplicative factor
    nLambda       = ceil(log(lamMax/lamMin)/log(filterSpacing))+1
    lambdaR       = convert(Array{Float64,1}, exp.(log(lamMin):log(filterSpacing):log(lamMin)+log(filterSpacing)*(nLambda-1)))
    nExtra        = nLambda - length(Filter_base)
    rR            = rMax*exp.(-collect(0:nExtra)*log(filterSpacing))# ranges corresponding to lambdaR plus nExtra 1 lags in convolution:
    #define kr domain fields here for lagged convolution
    ErK   = 1im*zeros(2*length(freqs),length(lambdaR)) #J0 in row 1 J1 in row 2
    EbK   = 1im*zeros(2*length(freqs),length(lambdaR)) #J0 in row 1 J1 in row 2
    HrK   = 1im*zeros(2*length(freqs),length(lambdaR)) #J0 in row 1 J1 in row 2
    HbK   = 1im*zeros(2*length(freqs),length(lambdaR)) #J0 in row 1 J1 in row 2
    EzK   = 1im*zeros(length(freqs),length(lambdaR)) #only J1
    HzK   = 1im*zeros(length(freqs),length(lambdaR)) #only J1
    #needed for spline computations in space domain
    ErBase = 1im*zeros(length(rR),length(freqs))
    EbBase = 1im*zeros(length(rR),length(freqs))
    HrBase = 1im*zeros(length(rR),length(freqs))
    HbBase = 1im*zeros(length(rR),length(freqs))
    EzBase = 1im*zeros(length(rR),length(freqs))
    HzBase = 1im*zeros(length(rR),length(freqs))

    for ikrho in eachindex(lambdaR), iFreq in eachindex(freqs)
         ErK[2*iFreq-1:2*iFreq,ikrho], EbK[2*iFreq-1:2*iFreq,ikrho],
         HrK[2*iFreq-1:2*iFreq,ikrho], HbK[2*iFreq-1:2*iFreq,ikrho],
         EzK[iFreq,ikrho], HzK[iFreq,ikrho] = getCSEM1DKernelsAnisoHED(lambdaR[ikrho],freqs[iFreq],zRx[1],zTx,z,rho)
    end #first time loop over freq

    for iFreq in eachindex(freqs)
    # Extract results for all rR ranges:
        for irR  in eachindex(rR)
              iShift = (1:length(Filter_base)) .+ irR .-1
              ErBase[irR,iFreq] = dot(ErK[2*iFreq-1,iShift], Filter_J0)/rR[irR] + dot(ErK[2*iFreq,iShift], Filter_J1)/rR[irR]^2

              EbBase[irR,iFreq] = dot(EbK[2*iFreq-1,iShift], Filter_J0)/rR[irR] + dot(EbK[2*iFreq,iShift], Filter_J1)/rR[irR]^2

              HrBase[irR,iFreq] = dot(HrK[2*iFreq-1,iShift], Filter_J0)/rR[irR] + dot(HrK[2*iFreq,iShift], Filter_J1)/rR[irR]^2

              HbBase[irR,iFreq] = dot(HbK[2*iFreq-1,iShift], Filter_J0)/rR[irR] + dot(HbK[2*iFreq,iShift], Filter_J1)/rR[irR]^2

              EzBase[irR,iFreq] = dot(EzK[iFreq,iShift], Filter_J1)/rR[irR]

              HzBase[irR,iFreq] = dot(HzK[iFreq,iShift], Filter_J1)/rR[irR]
        end

        #fit spline at computed receiver locations
        #Dierckx needs x to be increasing so end:-1:1

        #Er
        splReal = Spline1D(log10.(rR[end:-1:1]), real(ErBase[end:-1:1,iFreq]))
        splImag = Spline1D(log10.(rR[end:-1:1]), imag(ErBase[end:-1:1,iFreq]))
        #evaluate spline at required rRx
        Er[:,iFreq] = evaluate(splReal, log10.(rRx)) - im*evaluate(splImag, log10.(rRx))

        #Eb
        splReal = Spline1D(log10.(rR[end:-1:1]), real(EbBase[end:-1:1,iFreq]))
        splImag = Spline1D(log10.(rR[end:-1:1]), imag(EbBase[end:-1:1,iFreq]))
        #evaluate spline at required rRx
        Eb[:,iFreq] = evaluate(splReal, log10.(rRx)) - im*evaluate(splImag, log10.(rRx))

        #Hr
        splReal = Spline1D(log10.(rR[end:-1:1]), real(HrBase[end:-1:1,iFreq]))
        splImag = Spline1D(log10.(rR[end:-1:1]), imag(HrBase[end:-1:1,iFreq]))
        #evaluate spline at required rRx
        Hr[:,iFreq] = evaluate(splReal, log10.(rRx)) - im*evaluate(splImag, log10.(rRx))

        #Hb
        splReal = Spline1D(log10.(rR[end:-1:1]), real(HbBase[end:-1:1,iFreq]))
        splImag = Spline1D(log10.(rR[end:-1:1]), imag(HbBase[end:-1:1,iFreq]))
        #evaluate spline at required rRx
        Hb[:,iFreq] = evaluate(splReal, log10.(rRx)) - im*evaluate(splImag, log10.(rRx))

        #Ez
        splReal = Spline1D(log10.(rR[end:-1:1]), real(EzBase[end:-1:1,iFreq]))
        splImag = Spline1D(log10.(rR[end:-1:1]), imag(EzBase[end:-1:1,iFreq]))
        #evaluate spline at required rRx
        Ez[:,iFreq] = evaluate(splReal, log10.(rRx)) - im*evaluate(splImag, log10.(rRx))

        #Hz
        splReal = Spline1D(log10.(rR[end:-1:1]), real(HzBase[end:-1:1,iFreq]))
        splImag = Spline1D(log10.(rR[end:-1:1]), imag(HzBase[end:-1:1,iFreq]))
        #evaluate spline at required rRx
        Hz[:,iFreq] = evaluate(splReal, log10.(rRx)) - im*evaluate(splImag, log10.(rRx))

    end #loop over frequencies
    #done with lagged convolution

return Er, Eb, Hr, Hb, Ez, Hz
end #function with lagged convolution

function getCSEM1DanisoHED(freqs::Array{T,1}, rRx::Array{T,1}, zRx::T, zTx::T, z::Array{T,1}, rho::Array{T,2}, fnum::Int;
     RxAzim::Float64 = 0.0, TxDip::Float64 = 0.0) where T<:AbstractFloat
    Er = 1.0im*zeros(length(rRx),length(freqs)) #space domain
    #lagged convolution, good for all Rx depths same
    lambdaR = 1.
    rR = 1.
    rMax   = maximum(rRx)
    rMin   = minimum(rRx)
    lamMin = Filter_base[1]/rMax
    lamMax = Filter_base[end]/rMin
    filterSpacing = Filter_base[2] / Filter_base[1] # the filters are spaced by a multiplicative factor
    nLambda       = ceil(log(lamMax/lamMin)/log(filterSpacing))+1
    lambdaR       = convert(Array{Float64,1}, exp.(log(lamMin):log(filterSpacing):log(lamMin)+log(filterSpacing)*(nLambda-1)))
    nExtra        = nLambda - length(Filter_base)
    rR            = rMax*exp.(-collect(0:nExtra)*log(filterSpacing))# ranges corresponding to lambdaR plus nExtra 1 lags in convolution:
    #define kr domain fields here for lagged convolution
    ErK   = 1im*zeros(3*length(freqs),length(lambdaR)) #J0 in row 1 J1 in row 2, J1V from VED
    #needed for spline computations in space domain
    ErBase = 1im*zeros(length(rR),length(freqs))

    for ikrho in eachindex(lambdaR), iFreq in eachindex(freqs)
         ErK[3*iFreq-2:3*iFreq,ikrho], = getCSEM1DKernelsAnisoHED(lambdaR[ikrho],freqs[iFreq],zRx[1],zTx,z,rho)
    end #first time loop over freq

    # Rx Azimuths and Tx Dips
    cst = cosd(TxDip)  # needed for all HED terms
    sit = sind(TxDip)  # needed for all VED terms
    csb = cosd(RxAzim) # needs to be multiplied into all HED/VED terms inside the loop
    sib = sind(RxAzim) # can be multiplied into all non-VED terms outside the loop

    for iFreq in eachindex(freqs)
    # Extract results for all rR ranges:
        for irR  in eachindex(rR)
              iShift = (1:length(Filter_base)) .+ irR .-1
              ErBase[irR,iFreq] = cst*csb*dot(ErK[3*iFreq-2,iShift], Filter_J0)/rR[irR] +
                                  cst*csb*dot(ErK[3*iFreq-1,iShift], Filter_J1)/rR[irR]^2 +
                                  sit*    dot(ErK[3*iFreq,  iShift], Filter_J1)/rR[irR]
        end

        #fit spline at computed receiver locations
        #Dierckx needs x to be increasing so end:-1:1

        #Er
        splReal = Spline1D(log10.(rR[end:-1:1]), real(ErBase[end:-1:1,iFreq]))
        splImag = Spline1D(log10.(rR[end:-1:1]), imag(ErBase[end:-1:1,iFreq]))
        #evaluate spline at required rRx
        Er[:,iFreq] = evaluate(splReal, log10.(rRx)) - im*evaluate(splImag, log10.(rRx))

    end #loop over frequencies
    #done with lagged convolution
    return Er
end #function with lagged convolution

function allocateEr(freqs::Array{Float64,1}, rRx::Array{Float64,1})
    Er = 1.0im*zeros(length(rRx),length(freqs)) #space domain
    #lagged convolution, good for all Rx depths same
    lambdaR = 1.
    rR = 1.
    rMax   = maximum(rRx)
    rMin   = minimum(rRx)
    lamMin = Filter_base[1]/rMax
    lamMax = Filter_base[end]/rMin
    filterSpacing = Filter_base[2] / Filter_base[1] # the filters are spaced by a multiplicative factor
    nLambda       = ceil(log(lamMax/lamMin)/log(filterSpacing))+1
    lambdaR       = convert(Array{Float64,1}, exp.(log(lamMin):log(filterSpacing):log(lamMin)+log(filterSpacing)*(nLambda-1)))
    nExtra        = nLambda - length(Filter_base)
    rR            = rMax*exp.(-collect(0:nExtra)*log(filterSpacing))# ranges corresponding to lambdaR plus nExtra 1 lags in convolution:
    #define kr domain fields here for lagged convolution
    ErK   = 1im*zeros(2*length(freqs),length(lambdaR)) #J0 in row 1 J1 in row 2
    #needed for spline computations in space domain
    ErBase = 1im*zeros(length(rR),length(freqs))
    return Er, ErK, ErBase, lambdaR, rR
end

function useallocatedEr(Er::Array{Complex{Float64}, 2}, ErK::Array{Complex{Float64}, 2}, ErBase::Array{Complex{Float64}, 2},
                        lambdaR::Array{Float64, 1}, rR::Array{Float64, 1},
                        freqs::Array{Float64,1}, rRx::Array{Float64,1}, zRx::Float64,
                        zTx::Float64, z::Array{Float64,1}, rho::Array{Float64,2})

    for ikrho in eachindex(lambdaR), iFreq in eachindex(freqs)
         ErK[2*iFreq-1:2*iFreq,ikrho], = getCSEM1DKernelsAnisoHED(lambdaR[ikrho],freqs[iFreq],zRx[1],zTx,z,rho)
    end #first time loop over freq

    for iFreq in eachindex(freqs)
    # Extract results for all rR ranges:
        for irR  in eachindex(rR)
              iShift = (1:length(Filter_base)) .+ irR .-1
              ErBase[irR,iFreq] = dot(ErK[2*iFreq-1,iShift], Filter_J0)/rR[irR] + dot(ErK[2*iFreq,iShift], Filter_J1)/rR[irR]^2
        end
        #fit spline at computed receiver locations
        #Dierckx needs x to be increasing so end:-1:1
        #Er
        splReal = Spline1D(log10.(rR[end:-1:1]), real(ErBase[end:-1:1,iFreq]))
        splImag = Spline1D(log10.(rR[end:-1:1]), imag(ErBase[end:-1:1,iFreq]))
        #evaluate spline at required rRx
        Er[:,iFreq] = evaluate(splReal, log10.(rRx)) - im*evaluate(splImag, log10.(rRx))
    end #loop over frequencies
    #done with lagged convolution
end

function efficientEr(freqs::Array{Float64,1}, rRx::Array{Float64,1},
                    zRx::Float64, zTx::Float64, z::Array{Float64,1}, rho::Array{Float64,2})
    Er, ErK, ErBase, lambdaR, rR = allocateEr(freqs, rRx)
    useallocatedEr(Er, ErK, ErBase, lambdaR, rR, freqs, rRx, zRx, zTx, z, rho)
    return Er
end

function getCSEM1DanisoHED(freqs::Array{T,1}, rRx::Array{T,1}, zRx::Array{T,1}, zTx::T, z::Array{T,1}, rho::Array{T,2}) where T<:AbstractFloat
    #not tested simple FHT!
      Er = 1.0im*zeros(length(rRx),length(freqs)) #space domain
      Eb = 1.0im*zeros(length(rRx),length(freqs))
      Hr = 1.0im*zeros(length(rRx),length(freqs))
      Hb = 1.0im*zeros(length(rRx),length(freqs))
      Ez = 1.0im*zeros(length(rRx),length(freqs))
      Hz = 1.0im*zeros(length(rRx),length(freqs))

      #do normal fast Hankel transform with digital filters (slow)
      ErK = 1im*zeros(2,length(Filter_base))   #k_r domain, 1st row for J0 2nd for J1
      EbK = 1im*zeros(2,length(Filter_base))
      HrK = 1im*zeros(2,length(Filter_base))
      HbK = 1im*zeros(2,length(Filter_base))
      EzK = 1im*zeros(1,length(Filter_base)) #only J1
      HzK = 1im*zeros(1,length(Filter_base)) #only J1

      for iFreq in eachindex(freqs)
        # Loop over receivers
        for iRx in eachindex(rRx)
          # Get the lambda values for this particular range rRxEval:
          lambda =  Filter_base/rRx[iRx]
          # Get kernels for the particular lambda's
          for ikrho in eachindex(lambda)
            ErK[:,ikrho], EbK[:,ikrho],
            HrK[:,ikrho], HbK[:,ikrho],
            EzK[:,ikrho], HzK[:,ikrho] = getCSEM1DKernelsAnisoHED(lambda[ikrho],freqs[iFreq],zRx[iRx],zTx,z,rho)
          end
          # Don't forget to scale the J1 kernel by 1/r for Er,Eb,Hr,Hb
          ErK[2,:] = ErK[2,:]/rRx[iRx]
          EbK[2,:] = EbK[2,:]/rRx[iRx]
          HrK[2,:] = HrK[2,:]/rRx[iRx]
          HbK[2,:] = HbK[2,:]/rRx[iRx]
          # Multiply kernels by filter weights and sum them all up,
          # remembering to divide by the particular range:
          Er[iRx,iFreq] =  (( ErK[1,:]*Filter_J0 + ErK[2,:]*Filter_J1 )/rRx[iRx])[1]

          Eb[iRx,iFreq] =  (( EbK[1,:]*Filter_J0 + EbK[2,:]*Filter_J1 )/rRx[iRx])[1]

          Hr[iRx,iFreq] =  (( HrK[1,:]*Filter_J0 + HrK[2,:]*Filter_J1 )/rRx[iRx])[1]

          Hb[iRx,iFreq] =  (( HbK[1,:]*Filter_J0 + HbK[2,:]*Filter_J1 )/rRx[iRx])[1]

          Ez[iRx,iFreq] = (EzK[1,:]*Filter_J1/rRx[iRx])[1]

          Hz[iRx,iFreq] = (HzK[1,:]*Filter_J1/rRx[iRx])[1]

        end # loop over receivers
      end #loop over frequencies
      return Er, Eb, Hr, Hb, Ez, Hz
end #function with normal FHT, no lagged convolution

end #module
