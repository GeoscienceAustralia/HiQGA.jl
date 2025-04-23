using HiQGA, PyPlot, DelimitedFiles
includet("RDP.jl")
nlayers = 52
vmin, vmax = -2.5, 0.5
items = [item for item in walkdir("/home/547/ar0754/z67/ar0754/largeaem/production")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28351, 28351, 28351, 28350, 28350, 28351, 28351, 28350]
out = map(zip(src_dir, src_epsg)) do (dir, epsg)
    RDP.writegiantfilefromxyzrhodir(nlayers; src_dir=dir, src_epsg=epsg)
end
## write to file
fname = "/home/547/ar0754/z67/ar0754/largeaem/production/giantfile_phase_02.txt"
format = "lat[1] long[1] src_epsg[1] X[1] Y[1] Z[1] zboundaries[52] log10rholow[52] log10rhomid[52] log10rhohigh[52] phid_mean[1] phid_sdev[1]"
write(fname, format)
io = open(fname, "a")
writedlm(io, reduce(vcat, reduce(vcat, out)))
close(io)