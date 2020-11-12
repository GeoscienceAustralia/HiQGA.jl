module TEMPEST1DInversion
import ..AbstractOperator.get_misfit
using ..AbstractOperator, ..AEM_VMD_HMD, Statistics
using PyPlot, LinearAlgebra, ..CommonToAll, Random, DelimitedFiles

const μ₀ = 4*pi*1e-7

mutable struct Bfield<:Operator1D
    F          :: AEM_VMD_HMD.HField
    dataBx     :: Array{Float64, 1}
    dataBz     :: Array{Float64, 1}
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
    Rot_tx     :: Array{Float64,2}
	x_rx       :: Float64
	y_rx       :: Float64
end

function Bfield(;
				dataBx = zeros(0),
				dataBy = zeros(0),
				σx = zeros(0),
				σy = zeros(0),
				selectx = zeros(Bool, 0),
				selecty = zeros(Bool, 0),
				useML = false,
				nfixed = 1,
				times = [],
				ramp  = [],
				z = zeros(0),
				ρ = zeros(0),
				ndatax = 0,
				ndataz = 0,
				ntimesperdecade = 10
				nfreqsperdecade = 5
				nkᵣeval = 60
				zTx = -120
				zRx = -80
				xRx = 115.0
				yRx = 0,
				rx_roll = 0.,
			    rx_pitch = 0.,
			    rx_yaw = 0.,
			    tx_roll = 0.,
			    tx_pitch = 0.,
			    tx_yaw = 0.,
				order = "ypr")

@assert !isempty(times)
@assert(!isempty(ramp))
@assert zTx<0
@assert zRx>zTx # receiver below transmitter
@assert xRx<0 # receiver behind transmitter

Rot_rx = makerotationmatrix(order=order,yaw=tx_yaw, pitch=tx_pitch, roll=tx_roll)
Rot_tx = makerotationmatrix(order=order,yaw=tx_yaw, pitch=tx_pitch, roll=tx_roll)

F = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      ntimesperdecade = ntimesperdecade,
                      nfreqsperdecade = nfreqsperdecade,
					  nkᵣeval = nkᵣeval,
                      times  = times,
                      ramp   = ramp,
                      zTx    = zTx,
                      rRx    = sqrt(x_rx^2 + y_rx^2)
                      zRx    = zRx,
					  modelprimary = false,
					  getradialH = true,
					  getazimH = true,
					  provideddt = false)

Bfield(F, dataBx, dataBz, useML,σx, σz, z, nfixed, ρ, selectx, selectz,
		ndatax, ndataz, rx_roll, rx_pitch, rx_yaw, tx_roll, tx_pitch, tx_yaw,
		Rot_rx, Rot_tx, x_rx, y_rx)

end

function makerotationmatrix(;yaw=0.0,roll=0.0,pitch=0.0, order="lala", doinv = false)
    ordervector = ["ypr","yrp","rpy","ryp","pyr","pry"]
    @assert any(ordervector .== order) """not a valid order such as "ypr" """
    # All these matrices are transposes
    Rr = [ 1.           0            0
           0            cosd(roll)  -sind(roll)
           0            sind(roll)   cosd(roll)  ]

    Rp = [ cosd(pitch)  0           -sind(pitch)
           0            1            0.
           sind(pitch)  0            cosd(pitch) ]

    Ry = [ cosd(yaw)   -sind(yaw)    0
           sind(yaw)    cosd(yaw)    0
           0            0            1.          ]
    # multiplying above transposes from left to right
    # is applying them to a column vector from right to left
    # e.g., ypr is r(p(y(v)))
    if order     == "ypr"
        Rot = Ry*Rp*Rr
    elseif order == "yrp"
        Rot = Ry*Rr*Rp
    elseif order == "rpy"
        Rot = Rr*Rp*Ry
    elseif order == "ryp"
        Rot = Rr*Ry*Rp
    elseif order == "pyr"
        Rot = Rp*Ry*Rr
    else
        Rot = Rp*Rr*Ry
    end

    if doinv
        Rot'
    else
        Rot
    end
end

end
