using HiQGA, PyPlot, DelimitedFiles
using RDP
nlayers = 52
vmin, vmax = -2.5, 0.5
items = [item for item in walkdir("/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28353, 28354, 28351, 28354, 28354, 28354, 28354, 28352, 28352, 28352, 28352]
out = map(zip(src_dir, src_epsg)) do (dir, epsg)
    dst_dir = joinpath(dirname(dir),"GDA94_vtk")
    RDP.writegiantfilefromxyzrhodir(nlayers; src_dir=dir, dst_dir, src_epsg=epsg)
end
## write to file
fname = "/Users/anray/Documents/work/projects/largeaem/final_01/giantfile.txt"
format = "lat long src_epsg X Y Z zboundaries log10ρlow log10ρmid log10ρhigh ϕmean ϕsdev"
write(fname, format)
io = open(fname, "a")
writedlm(io, reduce(vcat, reduce(vcat, out)))
close(io)