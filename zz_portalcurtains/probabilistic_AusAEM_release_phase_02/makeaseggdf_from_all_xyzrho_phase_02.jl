using HiQGA, PyPlot
using RDP
nlayers = 52
rootdir = "/scratch/ns59/ar0754/probabilistic_AEM_phase_02_products/asseggdf"
## multiple surveys
items = [item for item in walkdir("/g/data/z67/ar0754/largeaem/production")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28351, 28351, 28351, 28350, 28350, 28351, 28351, 28350]
map(zip(src_dir, src_epsg)) do (sdir, epsg)
    @info "doing $sdir"
    parts = split(sdir,"/")
    fname = parts[end-2]*"_"*parts[end-1]
    isSkyTEM = any(contains.(parts, "SkyTEM"))
    nunames = isSkyTEM ? nothing : ["z_rx", "x_rx"]
    nuunits = isSkyTEM ? nothing : ["m", "m"]
    dst_dir = joinpath(rootdir, parts[end-2])
    RDP.writeaseggdffromxyzrho(nlayers; src_dir=sdir, dst_dir, fname, src_epsg=epsg, nunames, nuunits)
end
