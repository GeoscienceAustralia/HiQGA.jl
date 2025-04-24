using HiQGA, PyPlot
using RDP
src_epsg = 28354
nlayers = 52
src_dir = "/Users/anray/Documents/work/projects/curtainstuff/summaries"
dst_dir = "/Users/anray/Documents/work/projects/curtainstuff/ERC_01_vtk_GDA94"
RDP.writevtkfromxyzrhodir(nlayers; src_dir, dst_dir, src_epsg)