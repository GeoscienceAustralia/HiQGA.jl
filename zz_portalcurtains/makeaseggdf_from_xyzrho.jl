using HiQGA, PyPlot
includet("RDP.jl")
src_epsg = 28354
nlayers = 52
src_dir = "/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all/AusAEM_03_ERC/2021_ERC_01/summary"
dst_dir = "/Users/anray/Documents/work/projects/largeaem/final_01/aseggdf/"
parts = split(src_dir,"/")
fname = parts[end-2]*"_"*parts[end-1]
RDP.writeaseggdffromxyzrho(nlayers; src_dir, dst_dir, 
         fname, src_epsg)
