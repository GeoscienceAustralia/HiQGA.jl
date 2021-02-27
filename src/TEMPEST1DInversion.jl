module TEMPEST1DInversion
import ..AbstractOperator.get_misfit
using ..AbstractOperator, ..AEM_VMD_HMD, Statistics
using PyPlot, LinearAlgebra, ..CommonToAll, Random, DelimitedFiles

import ..Model, ..Options, ..OptionsStat, ..OptionsNonstat
import ..ModelNuisance, ..OptionsNuisance

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
				peakcurrent = 0.5)

	@assert !isempty(times)
	@assert(!isempty(ramp))
	@assert zTx<0
	if strictgeometry
		@assert zRx>zTx # receiver below transmitter
		@assert x_rx<0 # receiver behind transmitter
	end
	Rot_rx = makerotationmatrix(order=order_rx,yaw=rx_yaw, pitch=rx_pitch, roll=rx_roll,doinv=true)
	Rot_tx = makerotationmatrix(order=order_tx,yaw=tx_yaw, pitch=tx_pitch, roll=tx_roll)
	F = AEM_VMD_HMD.HFieldDHT(;
	                      ntimesperdecade = ntimesperdecade,
	                      nfreqsperdecade = nfreqsperdecade,
						  freqlow=1e-5,
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
	Bfield(F, dataHx, dataHz, useML,σx, σz, z, nfixed, copy(ρ), selectx, selectz,
			ndatax, ndataz, rx_roll, rx_pitch, rx_yaw, tx_roll, tx_pitch, tx_yaw,
			Rot_rx, x_rx, y_rx, mhat, Hx, Hy, Hz, addprimary, peakcurrent)
end

#TODO for nuisance moves in an MCMC chain
function update_geometry(tempest::Bfield, geovec::Array{Float64,1},
	order_rx = "ypr", order_tx = "ypr")
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

	#make new rotation matrices
	Rot_rx = makerotationmatrix(order = order_rx,
		yaw = rx_yaw, pitch = rx_pitch, roll = rx_roll, doinv = true)
	Rot_tx = makerotationmatrix(order = order_tx,
		yaw = tx_yaw, pitch = tx_pitch, roll = tx_roll)

	#do update on internal VMD model
	AEM_VMD_HMD.update_ZR!(tempest.F, zTx, zRx, nothing, sqrt(x_rx^2 + y_rx^2))
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

#match API for SkyTEM inversion getfield
function getfield!(m::Model, tempest::Bfield)
	getfield!(m.fstar, tempest)
	nothing
end
function getfield!(m::Array{Float64}, tempest::Bfield)
    copyto!(tempest.ρ, tempest.nfixed+1:length(tempest.ρ), 10 .^m, 1:length(m))
	getfieldTD!(tempest, tempest.z, tempest.ρ)
    nothing
end

# set the field given a conductivity model (GP parametrisation)
# and nuisance model (vector)
function getfield!(m::Model, mn::ModelNuisance, tempest::Bfield)
	update_geometry(tempest, mn.nuisance)
	getfield!(m, tempest)
end

function getfield!(m::Array{Float64}, mn::Array{Float64}, tempest::Bfield)
    copyto!(tempest.ρ, tempest.nfixed+1:length(tempest.ρ), 10 .^m, 1:length(m))
	update_geometry(tempest, vec(mn))
	getfieldTD!(tempest, tempest.z, tempest.ρ)
    nothing
end

function reducegreenstensor!(tempest)
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

	HMDx = [y2mx2/r3*J1h + x2/r2*J0v,
		    -2xy/r3*J1h  + xy/r2*J0v,
		    -x/r*J1v                       ]

	HMDy = [HMDx[2],
	 		-y2mx2/r3*J1h + y2/r2*J0v,
	 		-y/r*J1v		    		   ]

	VMDz = [x/r*J1v,
			y/r*J1v,
			J0v                            ]

	if tempest.addprimary
		HMDxp = [3x2 - R2, 3xy     , 3xz      ]/fpiR5
		HMDyp = [3xy     , 3y2 - R2, 3yz      ]/fpiR5
		VMDzp = [3xz     , 3yz     , 3z^2 - R2]/fpiR5
		for idim in 1:3
			HMDx[idim] .+= currentfac*HMDxp[idim]
			HMDy[idim] .+= currentfac*HMDyp[idim]
			VMDz[idim] .+= currentfac*VMDzp[idim]
		end
	end

	Hx[:], Hy[:], Hz[:] = Rot_rx*Roll180*[HMDx HMDy VMDz]*Roll180*mhat
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

function get_misfit(m::Model, opt::Options, tempest::Bfield)
	chi2by2 = 0.0;
	if !opt.debug
		getfield!(m, tempest)
		chi2by2 = getchi2by2([tempest.Hx; tempest.Hz], [tempest.dataHx; tempest.dataHz],
		[tempest.σx;tempest.σz], false, tempest.ndatax + tempest.ndataz);
	end
	return chi2by2
end

function get_misfit(m::Model, mn::ModelNuisance, opt::Union{Options,OptionsNuisance}, tempest::Bfield)
	chi2by2 = 0.0;
	if !opt.debug
		getfield!(m, mn, tempest)
		chi2by2 = getchi2by2([tempest.Hx; tempest.Hz], [tempest.dataHx; tempest.dataHz],
		[tempest.σx; tempest.σz], false, tempest.ndatax + tempest.ndataz);
	end
	return chi2by2
end


function getchi2by2(fm, d, σ, useML, ndata)
    r = (fm - d)./σ
    if useML
        chi2by2 = 0.5*ndata[ifreq]*log(norm(r[idx])^2)
    else
        chi2by2 = 0.5*norm(r)^2
    end
end

function plotmodelfield!(tempest::Bfield, z::Array{Float64,1}, ρ::Array{Float64,1})
	getfieldTD!(tempest, z, ρ)
	figure(figsize=(8,6))
	s = subplot(121)
	s.step(vcat(ρ[2:end],ρ[end]), vcat(z[2:end], z[end]+50))
	s.set_xscale("log")
	s.invert_xaxis()
	s.invert_yaxis()
	xlabel("ρ")
	ylabel("depth m")
	grid(true, which="both")
	s1 = subplot(122)
	semilogx(tempest.F.times, abs.(μ₀*tempest.Hz)*fTinv, label="Bz")
	semilogx(tempest.F.times, abs.(μ₀*tempest.Hx)*fTinv, label="Bx")
	if !isempty(tempest.σx)
		errorbar(tempest.F.times, μ₀*abs.(vec(tempest.dataHz))*fTinv, yerr = μ₀*2vec(tempest.σz)*fTinv,
                         linestyle="none", marker=".", elinewidth=0, capsize=3)
		errorbar(tempest.F.times, μ₀*abs.(vec(tempest.dataHx))*fTinv, yerr = μ₀*2vec(tempest.σx)*fTinv,
						 linestyle="none", marker=".", elinewidth=0, capsize=3)
	end
	xlabel("time s")
	ylabel("B field 10⁻¹⁵ T")
	legend()
	grid(true, which="both")
	!tempest.addprimary && s1.set_yscale("log")
	nicenup(gcf())
	nothing
end
function plotmodelfield!(tempest::Bfield, Ρ::Vector{Array{Float64}};
                        figsize=(8,5), dz=-1., onesigma=true,
                        extendfrac=-1., fsize=10, alpha=0.1)
	times, f, ax, nfixed, z  = setupaxis(tempest, Ρ, figsize, dz, onesigma,
                          extendfrac, fsize, alpha)
    for ρ in Ρ
        getfield!(ρ, tempest)
        tempest.Hz[.!tempest.selectz] .= NaN
        tempest.Hx[.!tempest.selectx] .= NaN
        ax[1].step(log10.(tempest.ρ[2:end]), tempest.z[2:end], "-k", alpha=alpha)
        ax[2].semilogx(times,μ₀*abs.(tempest.Hz)*fTinv, "k", alpha=alpha, markersize=2)
        ax[2].semilogx(times,μ₀*abs.(tempest.Hx)*fTinv, "k", alpha=alpha, markersize=2)
    end
	finishaxis(ax, f, z, dz, extendfrac, nfixed, fsize)
end

function plotmodelfield!(tempest::Bfield, Ρ::Vector{Array{Float64}},
						mn::Array{Float64,2};
                        figsize=(8,5), dz=-1., onesigma=true,
                        extendfrac=-1., fsize=10, alpha=0.1)
	@assert length(Ρ) == size(mn, 1)
	times, f, ax, nfixed, z  = setupaxis(tempest, Ρ, figsize, dz, onesigma,
                          extendfrac, fsize, alpha)
    for (i, ρ) in enumerate(Ρ)
        getfield!(ρ, mn[i,:], tempest)
        tempest.Hz[.!tempest.selectz] .= NaN
        tempest.Hx[.!tempest.selectx] .= NaN
        ax[1].step(log10.(tempest.ρ[2:end]), tempest.z[2:end], "-k", alpha=alpha)
        ax[2].semilogx(times,μ₀*abs.(tempest.Hz)*fTinv, "k", alpha=alpha, markersize=2)
        ax[2].semilogx(times,μ₀*abs.(tempest.Hx)*fTinv, "k", alpha=alpha, markersize=2)
    end
	finishaxis(ax, f, z, dz, extendfrac, nfixed, fsize)
end

function setupaxis(tempest::Bfield, Ρ,
                   figsize, dz, onesigma,
                   extendfrac, fsize, alpha)
    @assert all((dz, extendfrac) .> 0)
    sigma = onesigma ? 1.0 : 2.0
    f = figure(figsize=figsize)
    ax = Vector{PyPlot.PyObject}(undef, 3)
    ax[1] = subplot(121)
    ρmin, ρmax = extrema(vcat(Ρ...))
    delρ = ρmax - ρmin
    ax[1].set_xlim(ρmin-0.1delρ,ρmax+0.1delρ)
    nfixed, z = tempest.nfixed, tempest.z
    ax[1].plot([ρmin-0.1delρ,ρmax+0.1delρ], z[nfixed+1]*[1., 1], color="b")
    ax[2] = subplot(122)
	times = tempest.F.times
    Hz, Hx = tempest.Hz, tempest.Hx
    dataHx, σx = tempest.dataHx, tempest.σx
	dataHz, σz = tempest.dataHz, tempest.σz
    ax[2].errorbar(times, μ₀*abs.(dataHz)*fTinv, yerr = μ₀*sigma*σz*fTinv,
                        linestyle="none", marker=".", elinewidth=0, capsize=3, label="Bz")
    ax[2].errorbar(times, μ₀*abs.(dataHx)*fTinv, yerr = μ₀*sigma*σx*fTinv,
                        linestyle="none", marker=".", elinewidth=0, capsize=3, label="Bx")
	times, f, ax, nfixed, z
end

function finishaxis(ax, f, z, dz, extendfrac, nfixed, fsize)
	ax[1].grid()
    ax[1].set_ylabel("Depth m")
    ax[1].plot(xlim(), z[nfixed+1]*[1, 1], "--k")
    if dz > 0
        axn = ax[1].twinx()
        ax[1].get_shared_y_axes().join(ax[1],axn)
        yt = ax[1].get_yticks()[ax[1].get_yticks().>=z[nfixed+1]]
        axn.set_yticks(yt)
        axn.set_ylim(ax[1].get_ylim()[end:-1:1])
        axn.set_yticklabels(string.(Int.(round.(getn.(yt .- z[nfixed+1], dz, extendfrac)))))
    end
    axn.set_ylabel("Depth index", rotation=-90)
    ax[1].set_xlabel("Log₁₀ρ")
    ax[1].set_title("Model")
    ax[2].set_ylabel("B field 10¹⁵")
    ax[2].set_xlabel("Time (s)")
    ax[2].set_title("Transient response")
    ax[2].legend()
    ax[2].grid()
    ax[1].invert_xaxis()
    nicenup(f, fsize=fsize)
end

#for synthetics
function set_noisy_data!(tempest::Bfield, z::Array{Float64,1}, ρ::Array{Float64,1};
	noisefracx = 0.05, noisefracz = 0.05)
	primaryflag = tempest.addprimary
	if tempest.addprimary
		# adds noise only proportional to secondary field
		# when we compute full field next
		tempest.addprimary = false
	end
	getfieldTD!(tempest, z, ρ)
	σx = noisefracx*abs.(tempest.Hx)
	σz = noisefracz*abs.(tempest.Hz)
	# reset the tempest primary field modeling flag to original
	tempest.addprimary = primaryflag
	set_noisy_data!(tempest, z, ρ, σx, σz)
	plotmodelfield!(tempest, z, ρ)
	nothing
end

function set_noisy_data!(tempest::Bfield, z::Array{Float64,1}, ρ::Array{Float64,1},
	σx, σz;rseed = 123)
	Random.seed!(rseed)
	getfieldTD!(tempest, z, ρ)
	tempest.σx = σx
	tempest.σz = σz
	set_noisy_data(tempest,
		dataHx = tempest.Hx + σx.*randn(size(σx)),
		dataHz = tempest.Hz + σz.*randn(size(σz)),
		σx = σx,
		σz = σz)
	nothing
end

function set_noisy_data(tempest::Bfield;
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
	nothing
end

end
