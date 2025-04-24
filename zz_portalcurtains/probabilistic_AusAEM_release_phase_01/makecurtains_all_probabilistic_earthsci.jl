using HiQGA, PyPlot
using RDP
##
nlayers = 52
dr = 200
dz = 2
dpi = 400
fnamebar = "colorbar.jpg"
cmap = "turbo"
writegeom = true 
vmin, vmax = -2.5, 0.5
ϵfrac = 0.001
barfigsize = (0.4, 1.2)
shrink = 8000
VE = 30
rootdir = "/Users/anray/Documents/work/projects/curtainstuff/final_curtains_01/"
## multiple surveys
items = [item for item in walkdir("/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28353, 28354, 28351, 28354, 28354, 28354, 28354, 28352, 28352, 28352, 28352]
map(zip(src_dir, src_epsg)) do (sdir, epsg)
    prefix = split(sdir, "/")[end-2]*"_"*basename(dirname(sdir))*"_" 
    dst_dir = joinpath(rootdir, prefix[1:end-1]) # drop the ending underscore
    RDP.doallcurtaintriads.(;src_dir=sdir, dst_dir, nlayers, dr, dz, ϵfrac, src_epsg=epsg,
    barfigsize, dpi, cmap, fnamebar, shrink, VE, prefix,
    vmin, vmax, writegeom)
end    