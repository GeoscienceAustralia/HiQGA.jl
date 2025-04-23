using HiQGA, PyPlot
includet("RDP.jl")
nlayers = 52
rootdir = "/scratch/ns59/ar0754/probabilistic_AEM_phase_02_products/shapefiles/"
## multiple surveys
items = [item for item in walkdir("/g/data/z67/ar0754/largeaem/production")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28351, 28351, 28351, 28350, 28350, 28351, 28351, 28350]
map(zip(src_dir, src_epsg)) do (sdir, epsg)
    parts = split(sdir,"/")
    prefix = parts[end-2]*"_"*parts[end-1]
    RDP.writeshpfromxyzrhodir(nlayers; prefix, src_dir=sdir, dst_dir=rootdir, src_epsg=epsg)
end
