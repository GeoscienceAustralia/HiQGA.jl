using HiQGA, PyPlot
includet("RDP.jl")
##
nlayers = 52
dr = 480
dz = 2
dpi = 400
fnamebar = "colorbar.jpg"
cmap="turbo"
vmin, vmax = -0.5, 2.5
src_epsg = 28354
ϵfrac = 0.001
barfigsize = (0.4, 1.2)
shrink = 8000
VE = 30
prefix = "ERC_01"
## multiple lines
# destdir for multiple lines
dst_dir = "/Users/anray/Documents/work/projects/curtainstuff/ERC_01_curtains_for_mgn"
src_dir = "/Users/anray/Documents/work/projects/curtainstuff/summaries"
RDP.doallcurtaintriads.(;src_dir, dst_dir, nlayers, dr, dz, ϵfrac, src_epsg,
    barfigsize, dpi, cmap, fnamebar, shrink, VE, prefix,
    vmin, vmax)