using HiQGA, PyPlot
using RDP
nlayers = 52
rootdir = "/Users/anray/Documents/work/projects/largeaem/final_01/shapefiles/"
## multiple surveys
items = [item for item in walkdir("/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28353, 28354, 28351, 28354, 28354, 28354, 28354, 28352, 28352, 28352, 28352]
map(zip(src_dir, src_epsg)) do (sdir, epsg)
    parts = split(sdir,"/")
    prefix = parts[end-2]*"_"*parts[end-1]
    RDP.writeshpfromxyzrhodir(nlayers; prefix, src_dir=sdir, dst_dir=rootdir, src_epsg=epsg)
end
