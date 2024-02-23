using HiQGA, PyPlot
includet("RDP.jl")
nlayers = 52
items = [item for item in walkdir("/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28353, 28354, 28351, 28354, 28354, 28354, 28354, 28352, 28352, 28352, 28352]
map(zip(src_dir, src_epsg)) do (dir, epsg)
    dst_dir = joinpath(dirname(dir),"GDA94_vtk")
    RDP.writevtkfromxyzrhodir(nlayers; src_dir=dir, dst_dir, src_epsg=epsg)
end
## XML files
using HiQGA, Glob
idx_summary = [basename(it[1]) == "GDA94_vtk" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
vmin, vmax = -2.5, 0.5
map(src_dir) do (dir)
    fn = readdir(glob"Line*.vts", dir)
    # @info fn[1]
    # GDA94 to WGS84
    map(fn) do f
        transD_GP.writevtkxmlforcurtain(f; src_epsg=4283, dst_epsg=4326, suffix="", vmin, vmax)
    end    
end    