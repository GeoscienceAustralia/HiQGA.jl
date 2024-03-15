using Glob
rootdir = "/scratch/ns59/ar0754/"
## multiple surveys
items = [item for item in walkdir("/g/data/z67/ar0754/largeaem/production")]
idx_summary = [basename(it[1]) == "summary" for it in items]
src_dir = [it[1] for it in items[idx_summary]]
src_epsg = [28353, 28354, 28351, 28354, 28354, 28354, 28354, 28352, 28352, 28352, 28352]
rootdir = "/scratch/ns59/ar0754/probabilistic_AEM_phase_01_products"
map(zip(src_dir, src_epsg)) do (sdir, epsg)
    prefix = split(sdir, "/")[end-2]*"_"*basename(dirname(sdir))*"_" 
    containingdir = dirname(sdir)
    map(zip(["images", "sgrids", "vtk", "xyzlog10rho"],["pngs", "sgrids", "vtk", "summary"])) do (dstname, srcname)
        sd_use  = joinpath(containingdir, srcname)
        dst_dir = joinpath(rootdir, dstname, prefix[1:end-1])
        isdir(dst_dir) || mkpath(dst_dir)
        close(ioproj)
        if srcname == "summary"
            fn = readdir(glob"*_summary_xyzrho.txt", sd_use)
            cp.(fn, joinpath.(dst_dir, basename.(fn)))
        else
            cp(sd_use, dst_dir; force=true)
        end
        ioproj = open(joinpath(dst_dir, "0000_projection.txt"), "w")
        write(ioproj, "EPSG: $epsg")
        close(ioproj)
    end
    @info "done $sdir"
end