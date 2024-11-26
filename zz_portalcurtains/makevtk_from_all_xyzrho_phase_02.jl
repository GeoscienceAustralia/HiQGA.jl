using HiQGA, PyPlot
includet("RDP.jl")
nlayers = 52
vmin, vmax = -2.5, 0.5
outputrootdir = "/scratch/ns59/ar0754/probabilistic_AEM_phase_02_products"
items = [item for item in walkdir("/g/data/z67/ar0754/largeaem/production")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28351, 28351, 28351, 28350, 28350, 28351, 28351, 28350]
map(zip(src_dir, src_epsg)) do (dir, epsg)
    @show dst_dir = joinpath(outputrootdir, "GDA94_vtk", split(dir,"/")[end-1])
    RDP.writevtkfromxyzrhodir(nlayers; src_dir=dir, dst_dir, src_epsg=epsg, vmin, vmax)
end