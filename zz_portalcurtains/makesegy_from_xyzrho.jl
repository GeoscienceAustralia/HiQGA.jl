using HiQGA, PyPlot
includet("RDP.jl")
src_epsg = 28354
nlayers = 52
src_dir = "/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all/AusAEM_03_ERC/2021_ERC_01/summary"
dr = 250
dz = 1
parts = split(src_dir,"/")
fname = parts[end-2]*"_"*parts[end-1]
dst_dir = joinpath("/Users/anray/Documents/work/projects/largeaem/final_01/","segy",fname)
RDP.writesegyfromxyzrhodir(nlayers; src_dir, src_epsg, dst_dir, dr, dz)
