using HiQGA, PyPlot
cd(@__DIR__)
# This script writes out everything to GEODETIC GDA94
# other scripts for writing to SEG-Y or in local projections are in HiQGA.jl/zz_portalcurtains
# RDP is not a part of the HiQGA package but is an assemblage of functions put into a module on GitHub.
# Including it as part of HiQGA would mean dependencies on GDAL etc. which are best avoided
# to ensure HiQGA doesn't bloat
includet("/home/547/ar0754/.julia/dev/HiQGA/zz_portalcurtains/RDP.jl")
## probabilistic inversions
src_epsg = 28354
nlayers = 50
src_dir = "../UDF_probabilistic"
dst_dir = pwd()
# writes to dst_dir
RDP.writevtkfromxyzrhodir(nlayers; src_dir, dst_dir, src_epsg)
## deterministic inversions
fname = "../UDF_deterministic/UDF_deterministic_β²_0.1_R1_bg_0.01_Spm_zipped.dat"
epsg = 28354
dict = Dict("X"=>1, "Y"=>2, "Z"=>3, "cond"=>[62,111], "thick"=>[12,61], "line"=>5)
# writes out to the same directory as the ASEG-GDF file
RDP.colstovtk(dict, fname, epsg, decfactor=1, hasthick=false, islog10=true)
