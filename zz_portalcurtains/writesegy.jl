using SegyIO, PyPlot, HiQGA
cd(@__DIR__)
includet("RDP.jl")
pathname="/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all/AusAEM_03_ERC/2021_ERC_01/summary"
nlayers = 52
line = 1001001
X, Y, Z, zall, ρlow, ρmid, ρhigh, = transD_GP.readxyzrhoϕ(line, nlayers; pathname)
dr = 250.
dz = 2
RDP.XYZ_zmid_gridtoSEGY(-ρlow, X, Y, Z; dr, zall, dz)# one line

