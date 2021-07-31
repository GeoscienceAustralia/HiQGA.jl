module TEMPEST1DInversion
import ..AbstractOperator.get_misfit
import ..AbstractOperator.Sounding
import ..AbstractOperator.makeoperator
import ..AbstractOperator.make_tdgp_opt

using ..AbstractOperator, ..AEM_VMD_HMD, Statistics
using PyPlot, LinearAlgebra, ..CommonToAll, Random, DelimitedFiles

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
		idxx, idxz = tempest.selectx, tempest.selectz
        chi2by2 = getchi2by2([tempest.Hx[idxx]; tempest.Hz[idxz]],
                        [tempest.dataHx[idxx]; tempest.dataHz[idxz]],
                        [tempest.σx[idxx];tempest.σz[idxz]],
                        tempest.useML, tempest.ndatax + tempest.ndataz)
	end
	return chi2by2
end

function get_misfit(m::Model, mn::ModelNuisance, opt::Union{Options,OptionsNuisance}, tempest::Bfield)
	chi2by2 = 0.0;
	if !opt.debug
		getfield!(m, mn, tempest)
		idxx, idxz = tempest.selectx, tempest.selectz
        chi2by2 = getchi2by2([tempest.Hx[idxx]; tempest.Hz[idxz]],
                        [tempest.dataHx[idxx]; tempest.dataHz[idxz]],
                        [tempest.σx[idxx];tempest.σz[idxz]],
                        tempest.useML, tempest.ndatax + tempest.ndataz)
	end
	return chi2by2
end


function getchi2by2(fm, d, σ, useML, ndata)
    r = (fm - d)./σ
    if useML
        chi2by2 = 0.5*ndata*log(norm(r)^2)
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
	noisefracx = 0.02, noisefracz = 0.02,
	halt_X = nothing, halt_Z = nothing)
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
	if tempest.addprimary
		# adds noise only proportional to secondary field
		# when we compute full field next
		tempest.addprimary = false
	end
	getfieldTD!(tempest, z, ρ)
	σx = sqrt.((noisefracx*abs.(tempest.Hx)).^2 + (halt_X/μ₀).^2)
	σz = sqrt.((noisefracz*abs.(tempest.Hz)).^2 + (halt_Z/μ₀).^2)
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
	set_noisy_data!(tempest,
		dataHx = tempest.Hx + σx.*randn(size(σx)),
		dataHz = tempest.Hz + σz.*randn(size(σz)),
		σx = σx,
		σz = σz)
	nothing
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


function read_survey_files(;
    fname_dat="",
    fname_specs_halt="",
    frame_height = -2,
    frame_dz = -2,
    frame_dx = -2,
    frame_dy = -2,
	Hxp = -1,
	Hzp = -1,
    Hxs = [-2, 2],
    Hzs = [-2, 2],
    units = 1e-15,
	yaw_rx = -1,
	pitch_rx = -1,
	roll_rx = -1,
	yaw_tx = -1,
	pitch_tx = -1,
	roll_tx = -1,
    figsize = (16,11),
    makesounding = false,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.02,
    X = -1,
    Y = -1,
	Z = -1,
    fid = -1,
    linenum = -1,
	fsize = 13)

    @assert frame_height > 0
    @assert frame_dz > 0
    @assert frame_dx > 0
    @assert frame_dy > 0
    @assert all(Hxs .> 0)
    @assert all(Hzs .> 0)
	@assert Hxp > 0
	@assert Hzp > 0
    @assert X > 0
    @assert Y > 0
	@assert Z > 0
    @assert linenum > 0
    @assert fid > 0
	@assert pitch_rx > 0
	@assert roll_rx > 0
	@assert yaw_rx > 0
	@assert pitch_tx > 0
	@assert roll_tx > 0
	@assert yaw_tx > 0
	@assert 0 < multnoise < 1.0

    @info "reading $fname_dat"
    if dotillsounding!= nothing
        soundings = readdlm(fname_dat)[startfrom:skipevery:dotillsounding,:]
    else
        soundings = readdlm(fname_dat)[startfrom:skipevery:end,:]
    end
    easting = soundings[:,X]
    northing = soundings[:,Y]
	topo = soundings[:,Z]
    fiducial = soundings[:,fid]
    whichline = soundings[:,linenum]

	d_Hx = soundings[:,Hxs[1]:Hxs[2]] # secondary field
    d_Hz = soundings[:,Hzs[1]:Hzs[2]] # secondary field
    σ_Hx = multnoise*d_Hx # noise proportional to 2ndary
    σ_Hz = multnoise*d_Hz # noise proportional to 2ndary
	d_Hx .+= soundings[:,Hxp] # add primary field
	d_Hz .+= soundings[:,Hzp] # add primary field

    z_tx = soundings[:,frame_height]
    z_rx = -(z_tx + soundings[:,frame_dz]) # Flipping to my earth geometry
    z_tx = -z_tx # Flipping to my earth geometry
    x_rx = soundings[:,frame_dx]
	y_rx = -soundings[:,frame_dy] # Flipping to my earth geometry

	d_pitch_rx = -soundings[:,pitch_rx] # Flipping to GA-AEM geometry
	d_yaw_rx = -soundings[:,yaw_rx] # Flipping to GA-AEM geometry
	d_roll_rx = soundings[:,roll_rx]

	d_pitch_tx = -soundings[:,pitch_tx] # Flipping to GA-AEM geometry
	d_yaw_tx = -soundings[:,yaw_tx] # Flipping to GA-AEM geometry
	d_roll_tx = soundings[:,roll_tx]

    @info "reading $fname_specs_halt"
    include(fname_specs_halt)
    @assert size(d_Hx, 2) == length(times)
    @assert size(d_Hz, 2) == length(times)
    @assert size(d_Hx, 2) == length(Hx_add_noise)
    @assert size(d_Hz, 2) == length(Hz_add_noise)
	σ_Hx .*= units
	σ_Hz .*= units
    Hx_add_noise[:] .*= units
    Hz_add_noise[:] .*= units
    d_Hx[:]     .*= units
    d_Hz[:]     .*= -units # Flipping the Z component to align with GA_AEM rx
    f = figure(figsize=figsize)
    ax = Array{Any, 1}(undef, 4)
    ax[1] = subplot(2,2,1)
    nsoundings = size(soundings, 1)
    plot_dHx = permutedims(d_Hx)
    # plot_dHx[plot_dHx .<0] .= NaN
    pcolormesh(1:nsoundings, times, plot_dHx, shading="nearest")
    xlabel("sounding #")
    cbHx = colorbar()
    cbHx.ax.set_xlabel("Bx", fontsize=fsize)
    ylabel("time s")
    axHx = ax[1].twiny()
    axHx.semilogy(Hx_add_noise, times, "-k")
	axHx.semilogy(Hx_add_noise, times, "--w")
    axHx.set_xlabel("high alt noise")
    ax[2] = subplot(2,2,3,sharex=ax[1], sharey=ax[1])
    plot_dHz = permutedims(d_Hz)
    # plot_dHz[plot_dHM .<0] .= NaN
    pcolormesh(1:nsoundings, times, plot_dHz, shading="nearest")
    xlabel("sounding #")
    cbHz = colorbar()
    cbHz.ax.set_xlabel("Bz", fontsize=fsize)
    ylabel("time s")
    ax[2].invert_yaxis()
    axHz = ax[2].twiny()
    axHz.semilogy(Hz_add_noise, times, "-k")
	axHz.semilogy(Hz_add_noise, times, "--w")
    axHz.set_xlabel("high alt noise")
    ax[3] = subplot(2,2,2, sharex=ax[1])
	plot(1:nsoundings, z_tx, label="z_tx", "--")
	plot(1:nsoundings, z_rx, label="z_rx")
	plot(1:nsoundings, x_rx, label="x_rx")
	plot(1:nsoundings, y_rx, label="y_rx")
    legend()
    xlabel("sounding #")
    ylabel("earth frame geometry m")
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
    xlabel("sounding #")
    ylabel("rotations degrees")
	legend()
	nicenup(f, fsize=fsize)
    if makesounding
        s_array = Array{TempestSoundingData, 1}(undef, nsoundings)
        for is in 1:nsoundings
            l, fi = Int(whichline[is]), fiducial[is]
            @info "read $is out of $nsoundings"
            dHx, dHz = vec(d_Hx[is,:]), vec(d_Hz[is,:])
            σHx = sqrt.(vec(σ_Hx[is,:].^2) + Hx_add_noise.^2)
            σHz = sqrt.(vec(σ_Hz[is,:].^2) + Hz_add_noise.^2)
            s_array[is] = TempestSoundingData(
				"sounding_$(l)_$fi", easting[is], northing[is],
				topo[is], fi, l,
				x_rx[is], y_rx[is], z_rx[is],
                d_roll_rx[is], d_pitch_rx[is], d_yaw_rx[is],
				d_roll_tx[is], d_pitch_tx[is], d_yaw_tx[is],
				z_tx[is],
                times, ramp,
                σHx, σHz, dHx, dHz
                )
        end
        return s_array
    end
end

function makeoperator( sounding::TempestSoundingData;
                       zfixed   = [-1e5],
                       ρfixed   = [1e12],
                       zstart = 0.0,
                       extendfrac = 1.06,
                       useML = false,
                       dz = 2.,
                       ρbg = 10,
                       nlayers = 40,
                       ntimesperdecade = 10,
                       nfreqsperdecade = 5,
                       showgeomplot = false,
                       plotfield = false,
					   addprimary = true)
    @assert extendfrac > 1.0
    @assert dz > 0.0
    @assert ρbg > 0.0
    @assert nlayers > 1
    nmax = nlayers+1

    zall, znall, zboundaries = setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=showgeomplot)
    z, ρ, nfixed = makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
    ρ[z.>=zstart] .= ρbg
    ## Tempest operator creation from sounding data
	aem = Bfield(
    zTx = sounding.z_tx, zRx = sounding.z_rx,
	x_rx = sounding.x_rx, y_rx = sounding.y_rx,
    rx_roll = sounding.roll_rx, rx_pitch = sounding.pitch_rx, rx_yaw = sounding.yaw_rx,
    tx_roll = sounding.roll_tx, tx_pitch = sounding.pitch_tx, tx_yaw = sounding.yaw_tx,
	ramp = sounding.ramp, times = sounding.times, useML = useML,
	z=z, ρ=ρ,
	addprimary = addprimary #this ensures that the geometry update actually changes everything that needs to be
	)

	set_noisy_data!(aem,
		dataHx = sounding.Hx_data/μ₀,
		dataHz = sounding.Hz_data/μ₀,
		σx = sounding.σ_x/μ₀,
		σz = sounding.σ_z/μ₀)

	plotfield && plotmodelfield!(aem, z, ρ)

    aem, znall
end
function make_tdgp_opt(sounding::TempestSoundingData;
                    rseed = nothing,
                    znall = [1],
                    fileprefix = "sounding",
                    nmin = 2,
                    nmax = 40,
                    K = GP.Mat32(),
                    demean = false,
                    sampledc = true,
                    sddc = 0.01,
                    sdpos = 0.05,
                    sdprop = 0.05,
                    fbounds = [-0.5 2.5],
                    λ = [2],
                    δ = 0.1,
                    pnorm = 2,
                    save_freq = 50,
					nuisance_sdev   = [0.],
					nuisance_bounds = [0. 0.],
					updatenuisances = true,
                    C = nothing,
					dispstatstoscreen = false,
					restart = false
                    )
    sdev_dc = sddc*diff(fbounds, dims=2)[:]
    sdev_pos = [sdpos*abs(diff([extrema(znall)...])[1])]
    sdev_prop = sdprop*diff(fbounds, dims=2)[:]
    xall = permutedims(collect(znall))
    xbounds = permutedims([extrema(znall)...])

	history_mode = "w"
	restart && (history_mode = "a")

    updatenonstat = false
    needλ²fromlog = false
    if rseed != nothing
        Random.seed!(rseed)
    end

    opt = OptionsStat(fdataname = fileprefix*"_",
                            nmin = nmin,
                            nmax = nmax,
                            xbounds = xbounds,
                            fbounds = fbounds,
                            xall = xall,
                            λ = λ,
                            δ = δ,
                            demean = demean,
                            sampledc = sampledc,
                            sdev_dc = sdev_dc,
                            sdev_prop = sdev_prop,
                            sdev_pos = sdev_pos,
                            pnorm = pnorm,
                            quasimultid = false,
                            K = K,
                            save_freq = save_freq,
                            needλ²fromlog = needλ²fromlog,
                            updatenonstat = updatenonstat,
							updatenuisances = updatenuisances,
                            dispstatstoscreen = dispstatstoscreen,
							history_mode = history_mode
                            )

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

	optn = OptionsNuisance(opt;
                        sdev = copy(nuisance_sdev),
						bounds = bounds, C = C,
						updatenuisances = updatenuisances)

    opt, optn
end

function makeoperatorandoptions(soundings::Array{TempestSoundingData, 1};
                                                    nplot = 2,
                                                    zfixed   = [-1e5],
                                                    ρfixed   = [1e12],
                                                    zstart = 0.0,
                                                    extendfrac = 1.06,
                                                    useML = false,
                                                    dz = 2.,
                                                    ρbg = 10,
                                                    nlayers = 40,
                                                    ntimesperdecade = 10,
                                                    nfreqsperdecade = 5,
                                                    showgeomplot = false,
                                                    plotfield = true,
                                                    nmin = 2,
                                                    nmax = 40,
                                                    K = GP.Mat32(),
                                                    demean = false,
                                                    sampledc = true,
                                                    sddc = 0.01,
                                                    sdpos = 0.05,
                                                    sdprop = 0.05,
                                                    fbounds = [-0.5 2.5],
                                                    λ = [2],
                                                    δ = 0.1,
                                                    save_freq = 50,
                                                    nuisance_sdev   = [0.],
                                                    nuisance_bounds = [0. 0.],
                                                    C = nothing,
                                                    updatenuisances = true,
                                                    dispstatstoscreen = false
                                                    )
    for idx in randperm(length(soundings))[1:nplot]
            aem, znall = makeoperator(soundings[idx],
                                   zfixed = zfixed,
                                   ρfixed = ρfixed,
                                   zstart = zstart,
                                   extendfrac = extendfrac,
                                   dz = dz,
                                   ρbg = ρbg,
                                   nlayers = nlayers,
                                   ntimesperdecade = ntimesperdecade,
                                   nfreqsperdecade = nfreqsperdecade,
                                   showgeomplot = showgeomplot,
                                   plotfield = plotfield,
                                   useML = useML)
    
            opt, optn = make_tdgp_opt(soundings[idx],
                                    znall = znall,
                                    fileprefix = soundings[idx].sounding_string,
                                    nmin = nmin,
                                    nmax = nmax,
                                    K = K,
                                    demean = demean,
                                    sampledc = sampledc,
                                    sddc = sddc,
                                    sdpos = sdpos,
                                    sdprop = sdprop,
                                    fbounds = fbounds,
                                    save_freq = save_freq,
                                    λ = λ,
                                    δ = δ,
                                    nuisance_bounds = nuisance_bounds,
                                    nuisance_sdev = nuisance_sdev,
                                    updatenuisances = updatenuisances,
                                    C = C,
                                    dispstatstoscreen = dispstatstoscreen
                                    )
    end    

end

function summarypost(soundings::Array{TempestSoundingData, 1};
        qp1=0.05,
        qp2=0.95,
        burninfrac=0.5,
        zstart = 0.0,
        extendfrac = 1.06,
        useML = false,
        dz = 2.,
        nlayers = 40,
        nmin = 2,
        nmax = 40,
        K = GP.Mat32(),
        demean = false,
        sampledc = true,
        sddc = 0.01,
        sdpos = 0.05,
        sdprop = 0.05,
        fbounds = [-0.5 2.5],
        λ = [2],
        δ = 0.1,
        save_freq = 50,
        nuisance_sdev   = [0.],
        nuisance_bounds = [0. 0.],
        updatenuisances = true)

    @assert extendfrac > 1.0
    @assert dz > 0.0
    @assert nlayers > 1
    zall, znall, = setupz(zstart, extendfrac, dz=dz, n=nlayers)

    linename = "_line_$(soundings[1].linenum)_summary.txt"
    fnames = ["rho_low", "rho_mid", "rho_hi", "rho_avg",
        "ddz_mean", "ddz_sdev", "phid_mean", "phid_sdev",
        "nu_low", "nu_mid", "nu_high"].*linename

    idxnotzero = findidxnotzero(length(nuisance_sdev), nuisance_bounds)
    nunominal = zeros(length(idxnotzero), length(soundings))
    if isfile(fnames[1])
        @warn fnames[1]*" exists, reading stored values"
        pl, pm, ph, ρmean,
        vdmean, vddev, χ²mean, χ²sd, 
        nulow, numid, nuhigh = map(x->readdlm(x), fnames)
        for idx = 1:length(soundings)
            nunominal[:,idx] = getnufromsounding(soundings[idx])[idxnotzero]
        end                                 
    else
        # this is a dummy operator for plotting
        aem, = makeoperator(soundings[1])
        pl, pm, ph, ρmean, vdmean, vddev = map(x->zeros(length(zall), length(soundings)), 1:6)
        χ²mean, χ²sd = zeros(length(soundings)), zeros(length(soundings))
        nulow, numid, nuhigh  = map(x->zeros(length(idxnotzero), length(soundings)), 1:3)
        for idx = 1:length(soundings)
            opt, optn = make_tdgp_opt(soundings[idx],
                                        znall = znall,
                                        fileprefix = soundings[idx].sounding_string,
                                        nmin = nmin,
                                        nmax = nmax,
                                        K = K,
                                        demean = demean,
                                        sampledc = sampledc,
                                        sddc = sddc,
                                        sdpos = sdpos,
                                        sdprop = sdprop,
                                        fbounds = fbounds,
                                        save_freq = save_freq,
                                        λ = λ,
                                        δ = δ,
                                        nuisance_bounds = nuisance_bounds,
                                        nuisance_sdev = nuisance_sdev,
                                        updatenuisances = updatenuisances,
                                        dispstatstoscreen = false)
            opt.xall[:] .= zall
            @info "$idx out of $(length(soundings))\n"
            opt.xall[:] .= zall
            nunominal[:,idx] = mean(optn.bounds[idxnotzero,:], dims=2)
            pl[:,idx], pm[:,idx], ph[:,idx], ρmean[:,idx],
            vdmean[:,idx], vddev[:,idx] = CommonToAll.plot_posterior(aem, opt, burninfrac=burninfrac,
                                                    qp1=qp1, qp2=qp2,
                                                    doplot=false)
            h, nuquants = plot_posterior(aem, optn, burninfrac=burninfrac, doplot=false)
            nulow[:, idx]  .= nuquants[:, 1]
            numid[:, idx]  .= nuquants[:, 2]
            nuhigh[:, idx] .= nuquants[:, 3]
            χ² = 2*CommonToAll.assembleTat1(opt, :U, temperaturenum=1, burninfrac=burninfrac)
            ndata = sum(.!isnan.(soundings[idx].Hx_data)) +
                    sum(.!isnan.(soundings[idx].Hz_data))
            χ²mean[idx] = mean(χ²)/ndata
            χ²sd[idx]   = std(χ²)/ndata
            if useML 
                χ²mean[idx] -= log(ndata)
                χ²sd[idx]   -= log(ndata) # I think, need to check
            end    
        end
        # write in grid format
        for (fname, vals) in Dict(zip(fnames, [pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd, nulow, numid, nuhigh]))
            writedlm(fname, vals)
        end
        # write in x, y, z, rho format
        for (i, d) in enumerate([pl, pm, ph, ρmean])
            xyzrho = makearray(soundings, d, zall)
            writedlm(fnames[i][1:end-4]*"_xyzrho.txt", xyzrho)
        end
    end
    pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd, zall, nulow, numid, nuhigh, nunominal
end

function plotsummarygrids1(meangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname; qp1=0.05, qp2=0.95,
    figsize=(10,10), fontsize=12, cmap="viridis", vmin=-2, vmax=0.5, Eislast=true, Nislast=true,
    topowidth=2, idx=nothing, omitconvergence=false, useML=false, preferEright=false, preferNright=false,
    saveplot=false)

    f = figure(figsize=figsize)
    dr = diff(gridx)[1]
    f.suptitle(lname*" Δx=$dr m, Fids: $(length(R))", fontsize=fontsize)
    nrows = omitconvergence ? 4 : 5
    icol = 1
    s = Array{Any, 1}(undef, nrows)
    if !omitconvergence
        s[icol] = subplot(nrows, 1, icol)
        if useML
            plot(R, exp.(χ²mean))
            semilogy(R, ones(length(R)), "--k")
            ylabel("variance factor")
            title("Max likelihood variance adjustment")
        else
            plot(R, χ²mean)
            plot(R, ones(length(R)), "--k")
            fill_between(R, vec(χ²mean-χ²sd), vec(χ²mean+χ²sd), alpha=0.5)
            ylabel(L"ϕ_d")
            title("Data misfit")
        end
        icol += 1
    end
    s[icol] = omitconvergence ? subplot(nrows, 1, icol) : subplot(nrows, 1, icol, sharex=s[icol-1])
    imshow(plgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    title("Percentile $(round(Int, 100*qp1)) conductivity")
    ylabel("Height m")
    icol += 1
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1], sharey=s[icol-1])
    imshow(pmgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("Percentile 50 conductivity")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    ylabel("Height m")
    icol += 1
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1], sharey=s[icol-1])
    imshow(meangrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("Mean conductivity")
    ylabel("Height m")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    icol +=1
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1], sharey=s[icol-1])
    imlast = imshow(phgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    xlabel("Line distance m")
    title("Percentile $(round(Int, 100*qp2)) conductivity")
    ylabel("Height m")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    xlim(extrema(gridx))
    map(x->x.tick_params(labelbottom=false), s[1:end-1])
    map(x->x.grid(), s)
    nicenup(f, fsize=fontsize)
    plotNEWSlabels(Eislast, Nislast, gridx, gridz, s[(omitconvergence ? 1 : 2 ): end])
    f.subplots_adjust(bottom=0.125)
    cbar_ax = f.add_axes([0.125, 0.05, 0.75, 0.01])
    cb = f.colorbar(imlast, cax=cbar_ax, orientation="horizontal")
    cb.ax.set_xlabel("Log₁₀ S/m")
    (preferNright && !Nislast) && s[end].invert_xaxis()
    (preferEright && !Eislast) && s[end].invert_xaxis()
    saveplot && savefig(lname*".png", dpi=300)
end

function plotsummarygrids2(σmeangrid, ∇zmeangrid, ∇zsdgrid, cigrid, gridx, gridz, topofine, lname;
    qp1=0.05, qp2=0.95, Eislast=true, Nislast=true,
    figsize=(10,10), fontsize=12, cmap="viridis", vmin=-2, vmax=0.5, topowidth=2)

    f = figure(figsize=figsize)
    f.suptitle(lname, fontsize=fontsize)
    s1 = subplot(411)
    imshow(σmeangrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    title("Mean conductivity")
    ylabel("Height m")
    colorbar()
    s2 = subplot(412, sharex=s1)
    imshow(abs.(cigrid), cmap=cmap, aspect="auto",
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("CI: $(round(Int, 100*(qp2-qp1))) conductivity")
    ylabel("Height m")
    colorbar()
    s3 = subplot(413, sharex=s2, sharey=s2)
    imshow(∇zmeangrid, cmap=cmap, aspect="auto",
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("Mean conductivity vertical derivative")
    ylabel("Height m")
    colorbar()
    s4 = subplot(414, sharex=s2, sharey=s2)
    imlast = imshow(∇zsdgrid, cmap=cmap, aspect="auto",
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    xlabel("Line distance m")
    title("Std dev of conductivity vertical derivative")
    ylabel("Height m")
    colorbar()
    xlim(extrema(gridx))
    nicenup(f, fsize=fontsize)
end

function plotsummarygrids3(nuhigh, nulow, numid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname, nunominal, numsize=2; qp1=0.05, qp2=0.95,
    figsize=(10,10), fontsize=12, cmap="viridis", vmin=-2, vmax=0.5, Eislast=true, Nislast=true,
    topowidth=2, idx=nothing, useML=false, preferEright=false, preferNright=false, labelnu=[""],
    saveplot=false)
    
    # get the finely interpolated nuisances
    # nuhighfine, nulowfine, numidfine = map(X->gridpoints(R, gridx, X), [nuhigh, nulow, numid])

    f = figure(figsize=figsize)
    dr = diff(gridx)[1]
    f.suptitle(lname*" Δx=$dr m, Fids: $(length(R))", fontsize=fontsize)
    nrows = 4 + size(nulow, 1) # add the number of nuisances == no. of rows in nuhigh
    icol = 1
    s = Array{Any, 1}(undef, nrows)
    s[icol] = subplot(nrows, 1, icol)
    if useML
        plot(R, exp.(χ²mean))
        semilogy(R, ones(length(R)), "--k")
        ylabel("variance factor")
        title("Max likelihood noise variance adjustment")
    else
        plot(R, χ²mean)
        plot(R, ones(length(R)), "--k")
        fill_between(R, vec(χ²mean-χ²sd), vec(χ²mean+χ²sd), alpha=0.5)
        ylabel(L"ϕ_d")
        title("Data misfit")
    end
    icol += 1
    icol = plotnuquant(nulow, numid, nuhigh, nunominal, s, R, icol, nrows, numsize, labelnu)
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1])
    imshow(plgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    title("Percentile $(round(Int, 100*qp1)) conductivity")
    ylabel("Height m")
    icol += 1
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1], sharey=s[icol-1])
    imshow(pmgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    title("Percentile 50 conductivity")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    ylabel("Height m")
    icol += 1
    s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1], sharey=s[icol-1])
    imlast = imshow(phgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    plot(gridx, topofine, linewidth=topowidth, "-k")
    xlabel("Line distance m")
    title("Percentile $(round(Int, 100*qp2)) conductivity")
    ylabel("Height m")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    xlim(extrema(gridx))
    map(x->x.tick_params(labelbottom=false), s[1:end-1])
    map(x->x.grid(), s)
    nicenup(f, fsize=fontsize)
    plotNEWSlabels(Eislast, Nislast, gridx, gridz, s[size(nulow, 1)+1:end])
    f.subplots_adjust(bottom=0.125)
    cbar_ax = f.add_axes([0.125, 0.05, 0.75, 0.01])
    cb = f.colorbar(imlast, cax=cbar_ax, orientation="horizontal")
    cb.ax.set_xlabel("Log₁₀ S/m")
    (preferNright && !Nislast) && s[end].invert_xaxis()
    (preferEright && !Eislast) && s[end].invert_xaxis()
    saveplot && savefig(lname*".png", dpi=300)
end

function plotnuquant(nqlow, nqmid, nqhigh, nunominal, s, gridx, icol, nrows, ms=2, labelnu=[""])
    nnu = size(nqlow, 1)
    labelnu[1] == "" || @assert length(labelnu) == nnu
    for inu = 1:nnu
        s[icol] = subplot(nrows, 1, icol, sharex=s[icol-1])
        s[icol].fill_between(gridx, nqlow[inu,:], nqhigh[inu,:], alpha=0.5)
        s[icol].plot(gridx, nqmid[inu,:])
        s[icol].plot(gridx, nunominal[inu,:], "o", markersize=ms)
        labelnu[1] == "" || s[icol].set_title(labelnu[inu])
        icol += 1
    end
    icol    
end   

function summaryimages(soundings::Array{TempestSoundingData, 1};
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        zstart = 0.0,
                        extendfrac = 1.06,
                        useML = false,
                        dz = 2.,
                        nlayers = 40,
                        nmin = 2,
                        nmax = 40,
                        K = GP.Mat32(),
                        demean = false,
                        sampledc = true,
                        sddc = 0.01,
                        sdpos = 0.05,
                        sdprop = 0.05,
                        fbounds = [-0.5 2.5],
                        λ = [2],
                        δ = 0.1,
                        save_freq = 50,
                        nuisance_sdev   = [0.],
                        nuisance_bounds = [0. 0.],
                        updatenuisances = true,
                        dr = 10,
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap="viridis",
                        figsize=(6,10),
                        topowidth=2,
                        idx = nothing,
                        showderivs = false,
                        omitconvergence = false,
                        preferEright = false,
                        preferNright = false,
                        saveplot = false,
                        numsize=1,
                        labelnu = [""]
                        )
    @assert !(preferNright && preferEright) # can't prefer both labels to the right
    pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd, zall, 
    nulow, numid, nuhigh, nunominal = summarypost(soundings,
                                        qp1 = qp1,
                                        qp2 = qp2,
                                        burninfrac = burninfrac,
                                        zstart = zstart,
                                        extendfrac = extendfrac,
                                        useML = useML,
                                        dz = dz,
                                        nlayers = nlayers,
                                        nmin = nmin,
                                        nmax = nmax,
                                        K = K,
                                        demean = demean,
                                        sampledc = sampledc,
                                        sddc = sddc,
                                        sdpos = sdpos,
                                        sdprop = sdprop,
                                        fbounds = fbounds,
                                        λ = λ,
                                        δ = δ,
                                        save_freq = save_freq,
                                        nuisance_sdev   = nuisance_sdev,
                                        nuisance_bounds = nuisance_bounds,
                                        updatenuisances = updatenuisances)

    phgrid, plgrid, pmgrid, σmeangrid, ∇zmeangrid,
    ∇zsdgrid, gridx, gridz, topofine, R, Z = makesummarygrid(soundings, pl, pm, ph, ρmean,
                                            vdmean, vddev, zall, dz, dr=dr)

    lname = "Line $(soundings[1].linenum)"
    Eislast, Nislast = whichislast(soundings)
    plotsummarygrids1(σmeangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname, qp1=qp1, qp2=qp2,
        figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax, Eislast=Eislast,
        Nislast=Nislast, topowidth=topowidth, idx=idx, omitconvergence=omitconvergence, useML=useML,
        preferEright=preferEright, preferNright=preferNright, saveplot=saveplot)

    plotsummarygrids3(nuhigh, nulow, numid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname, nunominal, numsize, labelnu=labelnu, qp1=qp1, qp2=qp2,
    figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax, Eislast=Eislast,
    Nislast=Nislast, topowidth=topowidth, idx=idx, useML=useML,
    preferEright=preferEright, preferNright=preferNright, saveplot=saveplot)

    if showderivs
    cigrid = phgrid - plgrid
    plotsummarygrids2(σmeangrid, ∇zmeangrid, ∇zsdgrid, cigrid, gridx, gridz, topofine, lname, qp1=qp1, qp2=qp2,
        figsize=figsize, fontsize=fontsize, cmap=cmap, vmin=vmin, vmax=vmax, topowidth=topowidth,
        Eislast=Eislast, Nislast=Nislast)
    end
end

function plotindividualsoundings(soundings::Array{TempestSoundingData, 1};
                        burninfrac=0.5,
                        zstart = 0.0,
                        extendfrac = 1.06,
                        dz = 2.,
                        nlayers = 40,
                        nmin = 2,
                        nmax = 40,
                        K = GP.Mat32(),
                        demean = false,
                        sampledc = true,
                        sddc = 0.01,
                        sdpos = 0.05,
                        sdprop = 0.05,
                        fbounds = [-0.5 2.5],
                        λ = [2],
                        δ = 0.1,
                        save_freq = 50,
                        nuisance_sdev   = [0.],
                        nuisance_bounds = [0. 0.],
                        updatenuisances = true,
                        nbins=100,
                        figsize  = (12,6),
                        zfixed   = [-1e5],
                        ρfixed   = [1e12],
                        ntimesperdecade = 10,
                        nfreqsperdecade = 5,
                        computeforwards = false,
                        nforwards = 100,
                      idxcompute = [1])
    for idx = 1:length(soundings)
        if in(idx, idxcompute)
            @info "Sounding number: $idx"
            aem, znall = makeoperator(soundings[idx],
                zfixed = zfixed,
                ρfixed = ρfixed,
                zstart = zstart,
                extendfrac = extendfrac,
                dz = dz,
                nlayers = nlayers,
                ntimesperdecade = ntimesperdecade,
                nfreqsperdecade = nfreqsperdecade)
            opt, optn = make_tdgp_opt(soundings[idx],
                                        znall = znall,
                                        fileprefix = soundings[idx].sounding_string,
                                        nmin = nmin,
                                        nmax = nmax,
                                        K = K,
                                        demean = demean,
                                        sampledc = sampledc,
                                        sddc = sddc,
                                        sdpos = sdpos,
                                        sdprop = sdprop,
                                        fbounds = fbounds,
                                        save_freq = save_freq,
                                        λ = λ,
                                        δ = δ,
                                        nuisance_bounds = nuisance_bounds,
                                        nuisance_sdev = nuisance_sdev,
                                        updatenuisances = updatenuisances,
                                        dispstatstoscreen = false)
            zall, znall, = setupz(zstart, extendfrac, dz=dz, n=nlayers)    
            opt.xall[:] .= zall       
            getchi2forall(opt, alpha=0.8) # chi2 errors
            CommonToAll.getstats(opt) # ARs for GP model
            CommonToAll.getstats(optn) # ARs for nuisances
            plot_posterior(aem, opt, burninfrac=burninfrac, nbins=nbins, figsize=figsize) # GP models
            ax = gcf().axes
            ax[1].invert_xaxis()
            plot_posterior(aem, optn, burninfrac=burninfrac, nbins=nbins, figsize=figsize) # nuisances
            if computeforwards
                m = assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)
                mn = CommonToAll.assemblenuisancesatT(optn, temperaturenum=1, burninfrac=burninfrac)
                Random.seed!(10)
                plotmodelfield!(aem, m[randperm(length(m))[1:nforwards]],
                                                          mn[randperm(length(m))[1:nforwards],:],
                                                          dz=dz, extendfrac=extendfrac)
            end
        end
    end
end

end
