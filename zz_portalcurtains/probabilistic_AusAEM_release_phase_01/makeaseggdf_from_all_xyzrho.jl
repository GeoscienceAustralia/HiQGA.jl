using HiQGA, PyPlot
includet("RDP.jl")
nlayers = 52
rootdir = "/Users/anray/Documents/work/projects/largeaem/final_01/aseggdfwithnu/"
## multiple surveys
items = [item for item in walkdir("/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28353, 28354, 28351, 28354, 28354, 28354, 28354, 28352, 28352, 28352, 28352]
nunames = ["z_rx", "x_rx"]
nuunits = ["m", "m"]
map(zip(src_dir, src_epsg)) do (sdir, epsg)
    @info "doing $sdir"
    parts = split(sdir,"/")
    fname = parts[end-2]*"_"*parts[end-1]
    dst_dir = joinpath(rootdir, parts[end-2])
    RDP.writeaseggdffromxyzrho(nlayers; src_dir=sdir, dst_dir, fname, src_epsg=epsg, nunames, nuunits)
end
