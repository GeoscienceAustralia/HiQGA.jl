transD_GP.SurfaceRegression.plot_surface_posterior(opt, optlog10λ, cmapmean="RdBu", property="moho depth", 
            property_units="km", cmappdf="viridis")
## plot shapefile
ax = gcf().axes
using Shapefile, ArchGDAL, GeoInterface, GeoFormatTypes
t = Shapefile.Table("/home/547/ar0754/Downloads/STE_2016_AUST.shp")
for iaxis in 1:2
    for row in t
        if (row.STE_NAME16 == "Queensland") || (row.STE_NAME16 == "Northern Territory")
            xy = GeoInterface.coordinates(Shapefile.shape(row))
            for (i, points) in enumerate(xy)
                XY = permutedims(reverse(reduce(hcat, points[1]), dims=1))
                reprojXY = Vector{Array{Float64,1}}(undef, size(XY,1))
                for i in 1:size(XY, 1)
                        reprojXY[i] = XY[i,:]
                end
                XY = permutedims(reduce(hcat, ArchGDAL.reproject(reprojXY, EPSG(4283), EPSG(4462))))
                ax[iaxis].plot(XY[:,1], XY[:,2], "-k", linewidth=0.5)
            end
        end    
    end
end