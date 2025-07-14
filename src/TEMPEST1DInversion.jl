module TEMPEST1DInversion
import ..AbstractOperator.get_misfit, ..main
import ..AbstractOperator.Sounding
import ..AbstractOperator.makeoperator
import ..AbstractOperator.makebounds
import ..AbstractOperator.getoptnfromexisting
import ..AbstractOperator.getnufromsounding
import ..AbstractOperator.plotmodelfield!
import ..AbstractOperator.getresidual 
import ..AbstractOperator.returnforwrite
import ..AbstractOperator.getndata
import ..AbstractOperator.setnuboundsandstartforgradinv # for gradientbased
import ..AbstractOperator.setnuforinvtype

using ..AbstractOperator, ..AEM_VMD_HMD, ..SoundingDistributor
using PyPlot, LinearAlgebra, ..CommonToAll, Random, DelimitedFiles, Distributed, Dates, Statistics, SparseArrays

import ..Model, ..Options, ..OptionsStat, ..OptionsNonstat
import ..ModelNuisance, ..OptionsNuisance, ..findidxnotzero

const μ₀ = 4*pi*1e-7

mutable struct Bfield<:Operator1D
    F          :: AEM_VMD_HMD.HField
    dataHx     :: Array{Float64, 1}
    dataHz     :: Array{Float64, 1}
    useML      :: Bool
    σx         :: Array{Float64, 1}
    σz         :: Array{Float64, 1}
    z          :: Array{Float64, 1}
    nfixed     :: Int
    ρ          :: Array{Float64, 1}
    selectx    :: Array{Bool, 1}
    selectz    :: Array{Bool, 1}
    ndatax     :: Int
    ndataz     :: Int
    rx_roll    :: Float64
    rx_pitch   :: Float64
    rx_yaw     :: Float64
    tx_roll    :: Float64
    tx_pitch   :: Float64
    tx_yaw     :: Float64
    Rot_rx     :: Array{Float64,2}
	x_rx       :: Float64
	y_rx       :: Float64
	mhat       :: Array{Float64, 1}
	Hx         :: Array{Float64, 1}
	Hy         :: Array{Float64, 1}
	Hz         :: Array{Float64, 1}
	addprimary :: Bool
	peakcurrent:: Float64
    vectorsum  :: Bool
    J          :: AbstractArray
    W          :: SparseMatrixCSC
    res        :: Vector
end

# If needed to make z axis flip to align with GA-AEM
const Roll180 = [1. 0. 0.
				 0 -1  0.
				 0  0 -1]

const fTinv = 1e15

function Bfield(;
				dataHx = zeros(0),
				dataHz = zeros(0),
				σx = zeros(0),
				σz = zeros(0),
				selectx = zeros(Bool, 0),
				selectz = zeros(Bool, 0),
				useML = false,
				nfixed = 1,
				times = [],
				ramp  = [],
				z = zeros(0),
				ρ = zeros(0),
				ndatax = 0,
				ndataz = 0,
				ntimesperdecade = 10,
				nfreqsperdecade = 5,
				nkᵣeval = 60,
				zTx = -120,
				zRx = -80,
				x_rx = -115.0,
				y_rx = 0.,
				rx_roll = 0.,
			    rx_pitch = 0.,
			    rx_yaw = 0.,
			    tx_roll = 0.,
			    tx_pitch = 0.,
			    tx_yaw = 0.,
				order_tx = "ypr",
				order_rx = "ypr",
				strictgeometry = true,
				addprimary = false,
				peakcurrent = 0.5,
                vectorsum = false,
                calcjacobian = false)

	@assert !isempty(times)
	@assert(!isempty(ramp))
	@assert zTx<0
	if strictgeometry
		@assert zRx>zTx # receiver below transmitter
		@assert x_rx<0 # receiver behind transmitter
	end
    nmax = length(ρ)+1
	Rot_rx = makerotationmatrix(order=order_rx,yaw=rx_yaw, pitch=rx_pitch, roll=rx_roll,doinv=true)
	Rot_tx = makerotationmatrix(order=order_tx,yaw=tx_yaw, pitch=tx_pitch, roll=tx_roll)
	F = AEM_VMD_HMD.HFieldDHT(;nmax, calcjacobian,
	                      ntimesperdecade = ntimesperdecade,
	                      nfreqsperdecade = nfreqsperdecade,
						#   freqlow =1e-5,
						  nkᵣeval = nkᵣeval,
	                      times  = times,
	                      ramp   = ramp,
	                      zTx    = zTx,
	                      rRx    = sqrt(x_rx^2 + y_rx^2),
	                      zRx    = zRx,
						  modelprimary = false,
						  getradialH = true,
						  getazimH = true,
						  provideddt = false)
	mhat = Rot_tx*[0,0,1] # dirn cosines in inertial frame for VMDz
	Hx, Hy, Hz = map(x->zeros(size(times)), 1:3)
    # for Gauss-Newton
    res, J, W = allocateJ(length(times), nfixed, length(ρ), F.calcjacobian, vectorsum, σx, σz)
	Bfield(F, dataHx, dataHz, useML, σx, σz, z, nfixed, copy(ρ), selectx, selectz,
			ndatax, ndataz, rx_roll, rx_pitch, rx_yaw, tx_roll, tx_pitch, tx_yaw,
			Rot_rx, x_rx, y_rx, mhat, Hx, Hy, Hz, addprimary, peakcurrent, vectorsum, J, W, res)
end

function allocateJ(ntimes, nfixed, nmodel, calcjacobian, vectorsum, σx, σz)
    multfac = vectorsum ? 1 : 2
    if calcjacobian
        J = zeros(multfac*ntimes, nmodel-nfixed)
    else
        J = zeros(0)
    end
    # always return an allocated residuals and W - small price to pay I think
    # since majority of Jacobian allocations are in aem.F
    Wdiag = vectorsum ? ones(multfac*ntimes) : 1 ./[σx; σz]
    res = similar(Wdiag)
    W = sparse(diagm(Wdiag))
    res, J, W    
end

function allocateJ(tempest::Bfield)
    tempest.res, tempest.J, tempest.W = allocateJ(length(tempest.F.times), tempest.nfixed, 
        length(tempest.ρ),tempest.F.calcjacobian, tempest.vectorsum, tempest.σx, tempest.σz)
    nothing
end

function getresidual(tempest::Bfield, log10σ::Vector{Float64}, mn::Array{Float64}; computeJ=false) 
    # nuisance updates are separate to conductivity updates
    (; F, W, res) = tempest
    F.calcjacobian = computeJ
    getfield!(-log10σ, mn, tempest)
    # f = F.dBzdt[aem.select]
    # d = aem.d[aem.select]
    resforinvtype(tempest, W, res)
    nothing
end

function getresidual(tempest::Bfield, log10σ::Vector{Float64}; computeJ=false)
    (; F, W, res) = tempest
    F.calcjacobian = computeJ
    getfield!(-log10σ, tempest)
    # f = F.dBzdt[aem.select]
    # d = aem.d[aem.select]
    resforinvtype(tempest, W, res)
    nothing    
end    

function resforinvtype(tempest, W, res)
    if tempest.vectorsum
        f = get_fm(tempest)
        d, σ = get_dSigma(tempest)
        W[diagind(W)] = 1 ./σ
        res[:] = f - d
    else
        res[:] = [tempest.Hx; tempest.Hz] - [tempest.dataHx; tempest.dataHz]
    end 
end    

#TODO for nuisance moves in an MCMC chain
function update_geometry(tempest::Bfield, geovec::Array{Float64,1},
	order_rx = "ypr", order_tx = "ypr"; onlyprimary=false)
	length(geovec) == 10 ||
		throw(DimensionMismatch("TEMPEST geometry set with vector of wrong length."))
	zTx = geovec[1]
	zRx = geovec[2]
	x_rx = geovec[3]
	y_rx = geovec[4]

	rx_roll = geovec[5]
	rx_pitch = geovec[6]
	rx_yaw = geovec[7]

	tx_roll = geovec[8]
	tx_pitch = geovec[9]
	tx_yaw = geovec[10]

	# make new rotation matrices
	Rot_rx = makerotationmatrix(order = order_rx,
		yaw = rx_yaw, pitch = rx_pitch, roll = rx_roll, doinv = true)
	Rot_tx = makerotationmatrix(order = order_tx,
		yaw = tx_yaw, pitch = tx_pitch, roll = tx_roll)

	# do update on internal VMD model
    AEM_VMD_HMD.update_ZR!(tempest.F, zTx, zRx, nothing, sqrt(x_rx^2 + y_rx^2); onlyprimary)
	tempest.x_rx = x_rx
	tempest.y_rx = y_rx
	tempest.Rot_rx = Rot_rx
	tempest.mhat = Rot_tx*[0,0,1]

	tempest.rx_roll = rx_roll
	tempest.rx_pitch = rx_pitch
	tempest.rx_yaw = rx_yaw

	tempest.tx_roll = tx_roll
	tempest.tx_pitch = tx_pitch
	tempest.tx_yaw = tx_yaw

	nothing
end

function getfieldTD!(tempest::Bfield, z::Array{Float64, 1}, ρ::Array{Float64, 1})
	AEM_VMD_HMD.getfieldTD!(tempest.F,  z, ρ)
	reducegreenstensor!(tempest)
	nothing
end

function returnprimary!(tempest)
    # only useful for synthetics I guess
    geovec = extractnu(tempest)
    returnprimary!(tempest, geovec)
end

function returnprimary!(tempest, geovec)
    # only primary field calculations for given geovec
    update_geometry(tempest, geovec, onlyprimary=true)
    reducegreenstensor!(tempest, onlyprimary=true)
    tempest.Hx, tempest.Hy, tempest.Hz
end    

# all calling functions underneath here for misfit, field, etc. assume model is in log10 resistivity
# SANS the top. For lower level field calculation use AEM_VMD_HMD structs

function getfield!(m::Array{Float64}, tempest::Bfield)
    copyto!(tempest.ρ, tempest.nfixed+1:length(tempest.ρ), 10 .^m, 1:length(m))
	getfieldTD!(tempest, tempest.z, tempest.ρ)
    nothing
end

function getfield!(m::Array{Float64}, mn::Array{Float64}, tempest::Bfield)
    copyto!(tempest.ρ, tempest.nfixed+1:length(tempest.ρ), 10 .^m, 1:length(m))
	update_geometry(tempest, vec(mn))
	getfieldTD!(tempest, tempest.z, tempest.ρ)
    nothing
end

function reducegreenstensor!(tempest; onlyprimary = false)
	x, y   = tempest.x_rx, tempest.y_rx
	r      = tempest.F.rRx
	J1h    = tempest.F.dBazdt
	J1v    = tempest.F.dBrdt
	J0v    = tempest.F.dBzdt
	mhat   = tempest.mhat
	Rot_rx = tempest.Rot_rx
	Hx, Hy, Hz = tempest.Hx, tempest.Hy, tempest.Hz
	currentfac = tempest.peakcurrent

	y2mx2 = y^2-x^2
	r3    = r^3
	x2    = x^2
	r2    = r^2
	xy    = x*y
	y2    = y^2
	z     = tempest.F.zRx - tempest.F.zTx
	R2    = r2 + z^2
	fpiR5 = 4pi*sqrt(r2 + z^2)^5
	xz    = x*z
	yz    = y*z

    if !onlyprimary
        HMDx = [y2mx2/r3*J1h + x2/r2*J0v,
                -2xy/r3*J1h  + xy/r2*J0v,
                -x/r*J1v                       ]

        HMDy = [HMDx[2],
                -y2mx2/r3*J1h + y2/r2*J0v,
                -y/r*J1v		    		   ]

        VMDz = [x/r*J1v,
                y/r*J1v,
                J0v                            ]
    else
        HMDx, HMDy, VMDz = map(y->map(x->zeros(size(J0v)), 1:3), 1:3)
    end

    if tempest.addprimary || onlyprimary
        HMDxp = [3x2 - R2, 3xy     , 3xz      ]/fpiR5*currentfac
        HMDyp = [3xy     , 3y2 - R2, 3yz      ]/fpiR5*currentfac
        VMDzp = [3xz     , 3yz     , 3z^2 - R2]/fpiR5*currentfac
        for idim in 1:3
            if onlyprimary
                HMDx[idim] .= HMDxp[idim]
                HMDy[idim] .= HMDyp[idim]
                VMDz[idim] .= VMDzp[idim]
            else    
                HMDx[idim] .+= HMDxp[idim]
                HMDy[idim] .+= HMDyp[idim]
                VMDz[idim] .+= VMDzp[idim]
            end    
        end
    end

    Hx[:], Hy[:], Hz[:] = Rot_rx*Roll180*[HMDx HMDy VMDz]*Roll180*mhat        

    if tempest.F.calcjacobian
        F = tempest.F
        J_z, J_az, J_r = F.dBzdt_J, F.dBazdt_J, F.dBrdt_J
        HMDx_J = [y2mx2/r3*J_az + x2/r2*J_z,
                  -2xy/r3*J_az  + xy/r2*J_z,
                  -x/r*J_r                       ]

        HMDy_J = [HMDx_J[2],
                  -y2mx2/r3*J_az + y2/r2*J_z,
                  -y/r*J_r		    		     ]
   
        VMDz_J = [x/r*J_r,
                  y/r*J_r,
                  J_z                            ]
        
        Hx_J, _, Hz_J = Rot_rx*Roll180*[HMDx_J HMDy_J VMDz_J]*Roll180*mhat
        if tempest.vectorsum
            A = sqrt.(Hx.^2 + Hz.^2)        
            copy!(tempest.J, ((1 ./A').*((Hx').*Hx_J + (Hz').*Hz_J))'[:,tempest.nfixed+1:end-1])
        else
            copy!(tempest.J, hcat(Hx_J, Hz_J)'[:,tempest.nfixed+1:end-1])    
        end    
    end    
	nothing
end

function makerotationmatrix(;yaw=0.0,roll=0.0,pitch=0.0, order="lala", doinv = false)
    ordervector = ["ypr","yrp","rpy","ryp","pyr","pry"]
    @assert any(ordervector .== order) """not a valid order such as "pry" """
	# y --> z or [0,1,0] --> [0,0,1]
    Rr = [ 1.           0            0
           0            cosd(roll)  -sind(roll)
           0            sind(roll)   cosd(roll)  ]
	# z --> x or [0,0,1] --> [1,0,0]
    Rp = [ cosd(pitch)  0            sind(pitch)
           0            1            0.
          -sind(pitch)  0            cosd(pitch) ]
	# x --> y or [1,0,0] --> [0,1,0]
    Ry = [ cosd(yaw)   -sind(yaw)    0
           sind(yaw)    cosd(yaw)    0
           0            0            1.          ]
    # Applying them to a column vector from right to left
    # e.g., ypr is r(p(y(v)))
    if order     == "ypr"
        Rot = Rr*Rp*Ry
    elseif order == "yrp"
        Rot = Rp*Rr*Ry
    elseif order == "rpy"
        Rot = Ry*Rp*Rr
    elseif order == "ryp"
        Rot = Rp*Ry*Rr
    elseif order == "pyr"
        Rot = Rr*Ry*Rp
    else
        Rot = Ry*Rr*Rp
    end

    if doinv
        Rot'
    else
        Rot
    end
end

function get_misfit(m::AbstractArray, opt::Options, tempest::Bfield)
    # no nuisances
	chi2by2 = 0.0;
	if !opt.debug
		getfield!(m, tempest)
		chi2by2 = getchi2by2(tempest)
	end
	return chi2by2
end

function get_misfit(m::Model, opt::Options, tempest::Bfield)
    # actually used in McMC, no nuisances
    get_misfit(m.fstar, opt, tempest)
end

get_fm(tempest::Bfield) = get_fm(tempest.Hx, tempest.Hz)

function get_fm(Hx, Hz)
    fm = sqrt.(Hx.^2 + Hz.^2)
end 

get_dSigma(tempest::Bfield) = get_dSigma(tempest.dataHx,tempest.dataHz,tempest.σx,tempest.σz)

function get_dSigma(DHx, DHz, σx, σz)
    d = sqrt.(DHx.^2 + DHz.^2)
    σ = sqrt.((σx.^2).*DHx.^2 + (σz.^2).*DHz.^2)./d
    d, σ
end 

function get_misfit(m::AbstractArray, mn::AbstractVector, opt::Union{Options,OptionsNuisance}, tempest::Bfield)
    # this is useful when reading history files and debugging
	chi2by2 = 0.0;
	if !opt.debug
		getfield!(m, mn, tempest) # calls update_geometry down the stack
	    chi2by2 = getchi2by2(tempest)
	end
	return chi2by2
end

function get_misfit(σ::AbstractArray, nu::AbstractVector, tempest::Bfield, nubounds)
    # used by gradient inversion directly
    geovec = setnuforinvtype(tempest, nu) # mn is not as large as geovec
    lo, hi = nubounds[:,1], nubounds[:,2]
    (any(lo .> nu) | any(hi .< nu) ) && return 1e9
    tempest.F.calcjacobian = false
    getfield!(-σ, geovec, tempest) # calls update_geometry down the stack
    tempest.F.calcjacobian = true
	chi2by2 = getchi2by2(tempest)
end

function extractnu(tempest::Bfield)
    geovec = [
        tempest.F.zTx,
        tempest.F.zRx,
        tempest.x_rx,
        tempest.y_rx,
        tempest.rx_roll,
        tempest.rx_pitch,
        tempest.rx_yaw,
        tempest.tx_roll,
        tempest.tx_pitch,
        tempest.tx_yaw
    ]    
end    

function setnuforinvtype(tempest::Bfield, nu::AbstractVector)
    # order of elements must be right in nu!!
    geovec = extractnu(tempest)
    # required for amplitude only or vector sum
    geovec[2] = nu[1] # must be zRx
    geovec[3] = nu[2] # must be xRx
    if !tempest.vectorsum
        geovec[6] = nu[3] # receiver pitch
    end
    geovec
end    

function setnuboundsandstartforgradinv(tempest::Bfield, Δ::Array{Float64, 2})
    nubounds = similar(Δ)
    nustart = zeros(size(nubounds, 1))
    geovec = extractnu(tempest)
    # zRx and xRx low
    nubounds[1:2,1] = [geovec[2], geovec[3]] + Δ[1:2,1]
    nustart[1:2] = geovec[2:3]
    # zRx and xRx high
    nubounds[1:2,2] = [geovec[2], geovec[3]] + Δ[1:2,2]
    if !tempest.vectorsum
        nubounds[3,1] = geovec[6] + Δ[3,1]
        nubounds[3,2] = geovec[6] + Δ[3,2]
        nustart[3] = geovec[6]
    end
    nubounds, nustart    
end

function get_misfit(m::Model, mn::ModelNuisance, opt::Union{Options,OptionsNuisance}, tempest::Bfield)
    # actually used in McMC
    get_misfit(m.fstar, mn.nuisance, opt, tempest)
end

function getchi2by2(tempest::Bfield)
    idxx, idxz = tempest.selectx, tempest.selectz
    if tempest.vectorsum
        fm = get_fm(tempest)
        d, σ = get_dSigma(tempest)
        chi2by2 = getchi2by2(fm, d, σ, tempest.useML, tempest.ndatax)
    else
        chi2by2 = getchi2by2([tempest.Hx[idxx]; tempest.Hz[idxz]],
                    [tempest.dataHx[idxx]; tempest.dataHz[idxz]],
                    [tempest.σx[idxx];tempest.σz[idxz]],
                    tempest.useML, tempest.ndatax + tempest.ndataz)
    end
end

function getchi2by2(fm, d, σ, useML, ndata)
    r = (fm - d)./σ
    if useML
        chi2by2 = 0.5*ndata*log(norm(r)^2)
    else
        chi2by2 = 0.5*norm(r)^2
    end
end

# all plotting codes here assume that the model is in log10 resistivity, SANS
# the top layer resistivity. For lower level plotting use AEM_VMD_HMD structs

function plotdata(ax, d, σ, t; onesigma=true, dtype=nothing)
    sigma = onesigma ? 1 : 2
    color = "g"
    if dtype == :Hx
        label = "Bx"
        color = "m"
    elseif dtype == :Hz
        label = "Bz"
    else
        label = "|B|"
    end        
    ax.errorbar(t, μ₀*abs.(d)*fTinv; yerr = μ₀*sigma*fTinv*abs.(σ), color,
    linestyle="none", marker=".", elinewidth=1, capsize=3, label)
end

function plotsoundingcurve(ax, f, t; color="k", alpha=1, lw=1)
    ax.semilogx(t, μ₀*abs.(f)*fTinv, color=color, alpha=alpha, markersize=2, linewidth=lw)
end

function plotmodelfield!(ax, iaxis, aem::Bfield, ρ, nu; color=nothing, alpha=1, model_lw=1, forward_lw=1)
    # with nuisance
    stepmodel(ax, iaxis, color, ρ, aem, model_lw, alpha)
    getfield!(ρ, nu, aem)
    vectorsumsplit(ax, iaxis, aem, alpha, forward_lw, color)
end 

function plotmodelfield!(ax, iaxis, aem::Bfield, ρ; color=nothing, alpha=1, model_lw=1, forward_lw=1)
    # no nuisance
    stepmodel(ax, iaxis, color, ρ, aem, model_lw, alpha)
    getfield!(ρ, aem)
    vectorsumsplit(ax, iaxis, aem, alpha, forward_lw, color)
end 

function vectorsumsplit(ax, iaxis, aem::Bfield, alpha, forward_lw, color)
    colorused = !isnothing(color) ? color : ax[iaxis].lines[end].get_color()
    if aem.vectorsum
        fm = get_fm(aem)
        plotsoundingcurve(ax[iaxis+1], fm, aem.F.times; color=colorused, alpha, lw=forward_lw)
    else    
        plotsoundingcurve(ax[iaxis+1], aem.Hx, aem.F.times; color=colorused, alpha, lw=forward_lw)
        plotsoundingcurve(ax[iaxis+1], aem.Hz, aem.F.times; color=colorused, alpha, lw=forward_lw)
    end
end    

function initmodelfield!(aem;  onesigma=true, figsize=(8,6))
    f, ax = plt.subplots(1, 2; figsize)
    if !isempty(aem.dataHz)
        if (aem.vectorsum)
            d, σ = get_dSigma(aem)
            plotdata(ax[2], d, σ, aem.F.times; onesigma)
        else    
            aem.ndatax > 0 && plotdata(ax[2], aem.dataHx, aem.σx, aem.F.times; onesigma, dtype=:Hx)
            aem.ndataz > 0 && plotdata(ax[2], aem.dataHz, aem.σz, aem.F.times; onesigma, dtype=:Hz)
        end    
    end
    ax[1].set_xlabel("log10 ρ")
    ax[1].set_ylabel("depth m")
    ax[2].set_ylabel("B field fT")    
    ax[2].set_xlabel("time s")
    ax
end 

function plotmodelfield!(aem::Bfield, ρ::AbstractArray, nu::Vector; onesigma=true, color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,6), revax=true)
    # with nuisance
    plotmodelfield!(aem, [ρ], permutedims(nu); onesigma, color, alpha, model_lw, forward_lw, figsize, revax) 
end

function plotmodelfield!(aem::Bfield, ρ; onesigma=true, color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,6), revax=true)
    # no nuisance
    plotmodelfield!(aem, [ρ]; onesigma, color, alpha, model_lw, forward_lw, figsize, revax) 
end

function plotmodelfield!(aem::Bfield, manyρ::Vector{T}, manynu::Array{Float64, 2}; onesigma=true, 
        color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,6), revax=true) where T<:AbstractArray
    # with nuisance
    ax = initmodelfield!(aem; onesigma, figsize)
    for (i, ρ) in enumerate(manyρ)
        plotmodelfield!(ax, 1, aem, vec(ρ), manynu[i,:]; alpha, model_lw, forward_lw, color)
    end
    ax[1].invert_yaxis()
    nicenup(ax[1].get_figure(), fsize=12)
    revax && ax[1].invert_xaxis()
    ax
end

function plotmodelfield!(aem::Bfield, manyρ::Vector{T}; onesigma=true, 
        color=nothing, alpha=1, model_lw=1, forward_lw=1, figsize=(8,6), revax=true) where T<:AbstractArray
    # no nuisance    
    ax = initmodelfield!(aem; onesigma, figsize)
    for (i, ρ) in enumerate(manyρ)
        plotmodelfield!(ax, 1, aem, vec(ρ); alpha, model_lw, forward_lw, color)
    end
    ax[1].invert_yaxis()
    nicenup(ax[1].get_figure(), fsize=12)
    revax && ax[1].invert_xaxis()
    ax
end 

# for synthetics
function makenoisydata!(tempest::Bfield, ρ::Array{Float64,1};
	noisefracx = 0.02, noisefracz = 0.02, rseed=123, figsize=(8,5),
	halt_X = nothing, halt_Z = nothing, showplot=true,
    units = 1/fTinv)
	if halt_X != nothing
        @assert length(halt_X) == length(tempest.F.times)
    else
        halt_X = zeros(length(tempest.F.times))
    end
    if halt_Z != nothing
        @assert length(halt_Z) == length(tempest.F.times)
    else
        halt_Z = zeros(length(tempest.F.times))
    end
	primaryflag = tempest.addprimary
	# add noise only proportional to secondary field
    # when we compute full field next
	tempest.addprimary = false
	getfield!(ρ, tempest)
	σx = sqrt.((noisefracx*abs.(tempest.Hx)).^2 + (units*halt_X/μ₀).^2)
	σz = sqrt.((noisefracz*abs.(tempest.Hz)).^2 + (units*halt_Z/μ₀).^2)
	# reset the tempest primary field modeling flag to original
	tempest.addprimary = primaryflag
	# now compute full field with primary if flag says so (usual case)
    getfield!(ρ, tempest)
    Random.seed!(rseed)
    set_noisy_data!(tempest,
        dataHx = tempest.Hx + σx.*randn(size(σx)),
        dataHz = tempest.Hz + σz.*randn(size(σz)),
        σx = σx,
        σz = σz)
    
	if showplot
        plotwaveformgates(tempest)
        plotmodelfield!(tempest, ρ, figsize=figsize)
    end
    
	nothing
end

function makenoisydatafile!(fname::String, tempest::Bfield, ρ::Vector{Array{Float64,1}}, xrange;
	noisefracx = 0.02, noisefracz = 0.02,
	halt_X = nothing, halt_Z = nothing, 
    units = 1/fTinv)
    # remember to flip Hz, pitch and yaw from GA-AEM to z down TEMPEST
    # remember to flip y_rx and z_rx, z_tx from my system to TEMPEST
    d = map(enumerate(ρ)) do (i, rho)
        makenoisydata!(tempest, rho; 
            rseed=i, # clunky but ok
            noisefracx, noisefracz,
	        halt_X, halt_Z, showplot=false)
        Hxs, Hzs = copy(tempest.dataHx), copy(tempest.dataHz)     
        Hxp, _, Hzp = returnprimary!(tempest)
        Hxs, Hzs = Hxs-Hxp, Hzs-Hzp   
        hcat([i 1 xrange[i] 0 0 -tempest.F.zTx -abs(tempest.F.zTx-tempest.F.zRx)],
        [tempest.x_rx -tempest.y_rx],
        [tempest.rx_roll -tempest.rx_pitch -tempest.rx_yaw],
        [tempest.tx_roll -tempest.tx_pitch -tempest.tx_yaw],
        [Hxs'*μ₀/units -Hzs'*μ₀/units Hxp[1]*μ₀/units -Hzp[1]*μ₀/units])
    end
    reduce(vcat, d)
    headers = 
    """
    FID\t1
    Line\t2
    Easting\t3
    Northing\t4
    Height\t5
    frame_height\t6
    frame_dz\t7
    frame_dx\t8
    frame_dy\t9
    rx_roll\t10
    rx_pitch\t11
    rx_yaw\t12
    tx_roll\t13
    tx_pitch\t14
    tx_yaw\t15
    Hxs\t16-$(15+length(tempest.Hx))
    Hzs\t$(15+length(tempest.Hx)+1)-$(15+length(tempest.Hx)+length(tempest.Hz)) 
    Hxp\t$(15+length(tempest.Hx)+length(tempest.Hz)+1) 
    Hzp\t$(15+length(tempest.Hx)+length(tempest.Hz)+1+1)
    """
    f = open(fname*".hdr", "w")
    write(f, headers)
    close(f)
    writedlm(fname, d)
end

function set_noisy_data!(tempest::Bfield;
						dataHz = zeros(0), dataHx = zeros(0),
						σz = zeros(0), σx = zeros(0))

    @assert size(σx) == size(dataHx)
    @assert size(σz) == size(dataHz)
    ndatax  = sum(.!isnan.(dataHx))
    ndataz  = sum(.!isnan.(dataHz))
    selectx = .!isnan.(dataHx)
    selectz = .!isnan.(dataHz)

	tempest.σx = σx
	tempest.σz = σz
	tempest.dataHx = dataHx
	tempest.dataHz = dataHz
	tempest.selectx = selectx
	tempest.selectz = selectz
	tempest.ndatax = ndatax
	tempest.ndataz = ndataz

    # for Gauss-Newton
    allocateJ(tempest)

	nothing
end

mutable struct TempestSoundingData <: Sounding
    sounding_string :: String
    X       :: Float64
    Y       :: Float64
	Z       :: Float64
    fid     :: Float64
    linenum :: Int
    x_rx    :: Float64
	y_rx    :: Float64
	z_rx    :: Float64
	roll_rx :: Float64
	pitch_rx:: Float64
	yaw_rx  :: Float64
	roll_tx :: Float64
	pitch_tx:: Float64
	yaw_tx  :: Float64
    z_tx    :: Float64
    times   :: Array{Float64, 1}
    ramp    :: Array{Float64, 2}
    σ_x     :: Array{Float64, 1}
    σ_z     :: Array{Float64, 1}
    Hx_data :: Array{Float64, 1}
    Hz_data :: Array{Float64, 1}
    forceML :: Bool
    peakcurrent
    writefields 
end

function getnufromsounding(t::TempestSoundingData)
    # 1 zTx
    # 2 zRx
    # 3 x_rx
    # 4 y_rx
    # 5 rx_roll
    # 6 rx_pitch
    # 7 rx_yaw
    # 8 tx_roll
    # 9 tx_pitch
    # 10 tx_yaw
    return [t.z_tx, t.z_rx, t.x_rx, t.y_rx, t.roll_rx, t.pitch_rx, t.yaw_rx, 
            t.roll_tx, t.pitch_tx, t.yaw_tx]
end   

returnforwrite(s::TempestSoundingData) = getfield.(Ref(s), s.writefields)

function getndata(s::TempestSoundingData, vectorsum)
    n = getndata(s.Hx_data)[1] + getndata(s.Hz_data)[1]
    if vectorsum
        n = n/2 # tacitly assume number of X and Z are same 
    end
    n    
end    

function read_survey_files(;
    fname_dat="",
    fname_specs_halt="",
    frame_height = -99999999,
    frame_dz = -99999999,
    frame_dx = -99999999,
    frame_dy = -99999999,
	Hxp = -99999999,
	Hzp = -99999999,
    Hxs = [-99999999, -99999999],
    Hzs = [-99999999, -99999999],
    units = 1/fTinv,
	yaw_rx = -99999999,
	pitch_rx = -99999999,
	roll_rx = -99999999,
	yaw_tx = -99999999,
	pitch_tx = -99999999,
	roll_tx = -99999999,
    figsize = (10,6),
    dotillsounding = nothing,
    makeqcplots = true,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.02,
    peakcurrent = 0.5,
    X = -99999999,
    Y = -99999999,
	Z = -99999999,
    fid = -99999999,
    linenum = -99999999,
    lineslessthan = nothing,
	fsize = 10,
    dabsmin = nothing, # for very low amp
    earlyreturn = false, # to figure out very low amp,
    # these must now be passed in
    times = zeros(0),
    ramp = zeros(0,2),
    Hx_add_noise = zeros(0),
    Hz_add_noise = zeros(0),
    )

    @assert isempty(fname_specs_halt) """
        \n***
        include times, ramp, Hx_add_noise, Hz_add_noise explictly.
        If you already have an fname_specs_halt file, just include it
        before invoking read_survey_files() and use these variables
        as keyword args to this function
        ***
        """

    @assert frame_height > -99999999
    @assert frame_dz > -99999999
    @assert frame_dx > -99999999
    @assert frame_dy > -99999999
    @assert all(Hxs .> -99999999)
    @assert all(Hzs .> -99999999)
	@assert Hxp > -99999999
	@assert Hzp > -99999999
    @assert X > -99999999
    @assert Y > -99999999
	@assert Z > -99999999
    @assert linenum > -99999999
    @assert fid > -99999999
	@assert pitch_rx > -99999999
	@assert roll_rx > -99999999
	@assert yaw_rx > -99999999
	@assert pitch_tx > -99999999
	@assert roll_tx > -99999999
	@assert yaw_tx > -99999999
	@assert 0 < multnoise < 1.0
    if !isnothing(lineslessthan)
        lineslessthan::Int
    end    

    @info "reading $fname_dat"
    # soundings = readlargetextmatrix(fname_dat, startfrom, skipevery, dotillsounding)
    cols = [X, Y, Z, fid, linenum, [Hxs[1], Hxs[2]], [Hzs[1], Hzs[2]], Hxp, Hzp, frame_height,
    frame_dz, frame_dx, frame_dy, pitch_rx, yaw_rx, roll_rx, pitch_tx, yaw_tx, roll_tx]
    
    easting, northing, topo, fiducial, whichline, 
    d_Hxs, d_Hzs, d_Hxp, d_Hzp, z_tx, dz_rx, x_rx, y_rx, d_pitch_rx, d_yaw_rx, d_roll_rx,
    d_pitch_tx, d_yaw_tx, d_roll_tx = readcols(cols, fname_dat; decfactor=skipevery, 
                                        startfrom, dotill=dotillsounding)
    
    σ_Hx = multnoise*d_Hxs # noise proportional to 2ndary
    σ_Hz = multnoise*d_Hzs # noise proportional to 2ndary
	

    if frame_dz>0
        z_rx = -(z_tx + dz_rx) # Flipping to my earth geometry
        @warn "!!! flipping sign of z_rx to align with z down !!!"
    end
    if frame_height>0    
        z_tx = -z_tx # Flipping to my earth geometry
        @warn "!!! flipping sign of z_tx to align with z down !!!"
    end
    if frame_dy>0        
        y_rx = -y_rx # Flipping to my earth geometry
        @warn "!!! flipping sign of y_rx to align with z down !!!"
    end

    if pitch_rx>0
        d_pitch_rx = -d_pitch_rx # Flipping to GA-AEM geometry
        @warn "!!! flipping sign of pitch_rx to align with Rx z up !!!"
    end
    if yaw_rx>0    
        d_yaw_rx   = -d_yaw_rx # Flipping to GA-AEM geometry
        @warn "!!! flipping sign of yaw_rx to align with Rx z up !!!"
	end
    if pitch_tx>0
        d_pitch_tx = -d_pitch_tx # Flipping to GA-AEM geometry
        @warn "!!! flipping sign of pitch_tx to align with Rx z up !!!"
    end
    if yaw_tx>0    
        d_yaw_tx   = -d_yaw_tx # Flipping to GA-AEM geometry
        @warn "!!! flipping sign of yaw_tx to align with Rx z up !!!"
    end
    
    d_Hxs[:]     .*= units
    d_Hzs[:]     .*= units
    d_Hxp         *= units
    d_Hzp         *= units
    d_Hx = d_Hxs .+ d_Hxp
    @warn "!!! assuming Hzs and Hzp in same z dirn"
    d_Hz = d_Hzs .+ d_Hzp
    
    # this is now read explicitly
    # @info "reading $fname_specs_halt"
    # include(fname_specs_halt)
    @assert size(d_Hx, 2) == length(times)
    @assert size(d_Hz, 2) == length(times)
    @assert size(d_Hx, 2) == length(Hx_add_noise)
    @assert size(d_Hz, 2) == length(Hz_add_noise)
	σ_Hx .*= units
	σ_Hz .*= units
    Hx_add_noise[:] .*= units
    Hz_add_noise[:] .*= units
    if Hzp>0 || Hzs[1]>0
        d_Hz     = -d_Hz # Flipping the Z component to align with GA_AEM rx
        @warn "!!! flipping sign of d_Hz to align with Rx z up !!!"
    end    
    σ_Hx = sqrt.(σ_Hx.^2 .+ (Hx_add_noise').^2)
    σ_Hz = sqrt.(σ_Hz.^2 .+ (Hz_add_noise').^2)

    nsoundings = size(d_Hz, 1)
    makeqcplots && plotsoundingdata(nsoundings, times, d_Hxs, σ_Hx, d_Hzs, σ_Hz, z_tx, z_rx, x_rx, y_rx,
    d_yaw_tx, d_pitch_tx, d_roll_tx, d_yaw_rx, d_pitch_rx, d_roll_rx,
    figsize, fsize)

    earlyreturn && return d_Hxs, d_Hzs

    s_array = Array{TempestSoundingData, 1}(undef, nsoundings)
    fracdone = 0 
    countforceML = 0
    for is in 1:nsoundings
        l, fi = Int(whichline[is]), fiducial[is]
        if !isnothing(lineslessthan)
            l > lineslessthan && continue # skips high_alt and repeat lines if specified
        end    
        dHx, dHz = vec(d_Hx[is,:]), vec(d_Hz[is,:])
        forceML = checkifdatalow(d_Hxs[is,:], d_Hzs[is,:], dabsmin)
        countforceML += forceML ? 1 : 0
        s_array[is] = TempestSoundingData(
            "sounding_$(l)_$fi", easting[is], northing[is],
            topo[is], fi, l,
            x_rx[is], y_rx[is], z_rx[is],
            d_roll_rx[is], d_pitch_rx[is], d_yaw_rx[is],
            d_roll_tx[is], d_pitch_tx[is], d_yaw_tx[is],
            z_tx[is],
            times, ramp,
            σ_Hx[is,:], σ_Hz[is,:], dHx, dHz, forceML,
            peakcurrent, [:X, :Y, :Z, :fid, :linenum, :z_tx, :z_rx, :x_rx, :y_rx, :roll_rx, :pitch_rx, :yaw_rx, 
                :roll_tx, :pitch_tx, :yaw_tx]
            )
        fracnew = round(Int, is/nsoundings*100)
        if (fracnew-fracdone)>10
            fracdone = fracnew
            @info "read $is out of $nsoundings"
        end
    end
    idx = [isassigned(s_array, i) for i in 1:length(s_array)]
    !isnothing(dabsmin) && @info("low-amplitude forceML on $countforceML out of $(sum(idx)): $(round(100*countforceML/sum(idx)))%")
    return s_array[idx]
end

function checkifdatalow(d_Hxs, d_Hzs, dabsmin)
    forceML = false
    if !isnothing(dabsmin)
        dabs = sqrt.(d_Hxs.^2 + d_Hzs.^2)
        if mean(dabs) < dabsmin
            forceML = true
        end    
    end
    forceML    
end

function plotsoundingdata(nsoundings, times, d_Hx, Hx_add_noise, d_Hz, Hz_add_noise, z_tx, z_rx, x_rx, y_rx,
        d_yaw_tx, d_pitch_tx, d_roll_tx, d_yaw_rx, d_pitch_rx, d_roll_rx,
        figsize, fsize)
    # d_Hx, d_Hz are secondary fields!
    f = figure(figsize=figsize)
    ax = Array{Any, 1}(undef, 4)
    ax[1] = subplot(2,2,1)
    plot_dHx = permutedims(d_Hx)
    # plot_dHx[plot_dHx .<0] .= NaN
    im1 = ax[1].pcolormesh(1:nsoundings, times, plot_dHx, shading="nearest")
    ax[1].set_xlabel("sounding #")
    cbHx = colorbar(im1, ax=ax[1])
    cbHx.ax.set_xlabel("Bx secondary", fontsize=fsize)
    ax[1].set_ylabel("time s")
    axHx = ax[1].twiny()
    axHx.semilogy(infmean(Hx_add_noise./abs.(d_Hx), 1)[:], times, "r")
    axHx.semilogy(infmean(Hx_add_noise./abs.(d_Hx), 1)[:], times, "--w")
    axHx.set_xlabel("avg Hx noise fraction")
    ax[2] = subplot(2,2,3,sharex=ax[1], sharey=ax[1])
    plot_dHz = permutedims(d_Hz)
    # plot_dHz[plot_dHM .<0] .= NaN
    im2 = ax[2].pcolormesh(1:nsoundings, times, plot_dHz, shading="nearest")
    ax[2].set_xlabel("sounding #")
    cbHz = colorbar(im2, ax=ax[2])
    cbHz.ax.set_xlabel("Bz secondary", fontsize=fsize)
    ax[1].set_ylabel("time s")
    ax[2].set_ylabel("time s")
    ax[2].invert_yaxis()
    axHz = ax[2].twiny()
    axHz.semilogy(infmean(Hz_add_noise./abs.(d_Hz), 1)[:], times, "r")
    axHz.semilogy(infmean(Hz_add_noise./abs.(d_Hz), 1)[:], times, "--w")
    axHz.set_xlabel("avg Hz noise fraction")
    ax[3] = subplot(2,2,2, sharex=ax[1])
	ax[3].plot(1:nsoundings, z_tx, label="z_tx", "--")
	ax[3].plot(1:nsoundings, z_rx, label="z_rx")
	ax[3].plot(1:nsoundings, x_rx, label="x_rx")
	ax[3].plot(1:nsoundings, y_rx, label="y_rx")
    ax[3].legend()
    ax[3].set_xlabel("sounding #")
    ax[3].set_ylabel("earth frame geometry m")
	ax[3].grid()
    ax[3].invert_yaxis()
    ax[4] = subplot(2,2,4, sharex=ax[1])
    ax[4].plot(1:nsoundings, d_yaw_tx, label="yaw_tx", "--")
	ax[4].plot(1:nsoundings, d_pitch_tx, label="pitch_tx", "--")
	ax[4].plot(1:nsoundings, d_roll_tx, label="roll_tx", "--")
	ax[4].plot(1:nsoundings, d_yaw_rx, label="yaw_rx")
	ax[4].plot(1:nsoundings, d_pitch_rx, label="pitch_rx")
	ax[4].plot(1:nsoundings, d_roll_rx, label="roll_rx")
	ax[4].grid()
    ax[4].set_xlabel("sounding #")
    ax[4].set_ylabel("rotations degrees")
	ax[4].legend()
	nicenup(f, fsize=fsize)
end

function makeoperator( sounding::TempestSoundingData;
                       zfixed   = [-1e5],
                       ρfixed   = [1e12],
                       zstart = 0.0,
                       extendfrac = 1.06,
                       useML = false,
                       dz = 2.,
                       ρbg = 10,
                       calcjacobian = false, 
                       nlayers = 40,
                       ntimesperdecade = 10,
                       nfreqsperdecade = 5,
                       showgeomplot = false,
                       plotfield = false,
					   addprimary = true,
                       vectorsum = false)
    @assert extendfrac > 1.0
    @assert dz > 0.0
    @assert ρbg > 0.0
    @assert nlayers > 1

    zall, znall, zboundaries = setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=showgeomplot)
    z, ρ, = makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
    ρ[z.>=zstart] .= ρbg
    ## Tempest operator creation from sounding data
	aem = Bfield(;ntimesperdecade, nfreqsperdecade, calcjacobian,
    zTx = sounding.z_tx, zRx = sounding.z_rx,
	x_rx = sounding.x_rx, y_rx = sounding.y_rx,
    rx_roll = sounding.roll_rx, rx_pitch = sounding.pitch_rx, rx_yaw = sounding.yaw_rx,
    tx_roll = sounding.roll_tx, tx_pitch = sounding.pitch_tx, tx_yaw = sounding.yaw_tx,
	ramp = sounding.ramp, times = sounding.times, useML = useML,
	z=z, ρ=ρ,
	addprimary = addprimary, #this ensures that the geometry update actually changes everything that needs to be
    vectorsum = vectorsum,
    peakcurrent = sounding.peakcurrent
	)

	set_noisy_data!(aem,
		dataHx = sounding.Hx_data/μ₀,
		dataHz = sounding.Hz_data/μ₀,
		σx = sounding.σ_x/μ₀,
		σz = sounding.σ_z/μ₀)

	if plotfield 
        plotwaveformgates(aem)
        plotmodelfield!(aem, log10.(ρ[2:end]))
    end

    aem, zall, znall, zboundaries
end

function plotwaveformgates(aem::Bfield; figsize=(5,5))
    figure(;figsize)
    (;ramp, times) = aem.F
    plot(ramp[:,1]*1e6, ramp[:,2], "-or")
    stem(times*1e6, ones(length(times)))
    ylabel("Amplitude")
    xlabel("time μs")
    title("Ramp and gates linear time")
end  

function makeoperator(aem::Bfield, sounding::TempestSoundingData)
    ntimesperdecade = gettimesperdec(aem.F.interptimes)
    nfreqsperdecade = gettimesperdec(aem.F.freqs)
    @assert !(aem.useML & sounding.forceML) "useML and forceML through dabsmin cannot both be true"
    aemout = Bfield(;ntimesperdecade, nfreqsperdecade,
    zTx = sounding.z_tx, zRx = sounding.z_rx,
	x_rx = sounding.x_rx, y_rx = sounding.y_rx,
    rx_roll = sounding.roll_rx, rx_pitch = sounding.pitch_rx, rx_yaw = sounding.yaw_rx,
    tx_roll = sounding.roll_tx, tx_pitch = sounding.pitch_tx, tx_yaw = sounding.yaw_tx,
	ramp = sounding.ramp, times = sounding.times, 
    useML = (aem.useML | sounding.forceML), # OR logic for useML with true, true disqualified earlier
	z=copy(aem.z), ρ=copy(aem.ρ), calcjacobian=aem.F.calcjacobian,
	addprimary = aem.addprimary, #this ensures that the geometry update actually changes everything that needs to be
    vectorsum = aem.vectorsum,
    peakcurrent = sounding.peakcurrent
	)
    set_noisy_data!(aemout,
		dataHx = sounding.Hx_data/μ₀,
        dataHz = sounding.Hz_data/μ₀,
        σx = sounding.σ_x/μ₀,
        σz = sounding.σ_z/μ₀)
    aemout
end

function makebounds(nuisance_bounds, sounding::TempestSoundingData)
    bounds = nuisance_bounds .+ [
      sounding.z_tx
      sounding.z_rx
      sounding.x_rx
      sounding.y_rx
      sounding.roll_rx
      sounding.pitch_rx
      sounding.yaw_rx
      sounding.roll_tx
      sounding.pitch_tx
      sounding.yaw_tx]
end    

function getoptnfromexisting(optn_in::OptionsNuisance, opt_in::Options, sounding::TempestSoundingData)
    # helper function to get nuisances from an existing nuisance option
    (;sdev, idxnotzero, rotatebounds, bounds, idxnotzero) = optn_in
    sdevproposalfracs = zeros(size(sdev))
    sdevproposalfracs[idxnotzero] = sdev[idxnotzero]./diff(rotatebounds, dims=2)
    # assumption of symmetric priors about nuisance value!
    Δprior = 0.5diff(bounds, dims=2).*[-1 1]
    bounds = makebounds(Δprior, sounding)
    W = deepcopy(optn_in.W)
    OptionsNuisance(opt_in;
        sdev = copy(sdevproposalfracs),
        bounds, W,
        updatenuisances = true)
end  

end
