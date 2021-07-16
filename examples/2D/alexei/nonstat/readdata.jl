# this is to read alexei's files
using DelimitedFiles, PyPlot
flist = readdlm("filelist", String)
dvec = Array{Array{Float64}, 1}(undef, length(flist))
for (i,f) in enumerate(flist)
    dvec[i] = readdlm(f, Float64, skipstart=1)
end    
d = reduce(vcat, dvec)
long = d[:,1]
lat  = d[:,2]
σ = 1 ./d[:,3]
z = d[:,4]
# but I had to write it and reproject in Geosoft
writedlm("d_GDA94.txt", d)
# convert to Lambert conformal conic and then we have
d, h = readdlm("mohodepth_GDA_LCC.csv",',', header=true)
x, y = d[:,5], d[:,6]
f = figure()
s1 = subplot(121)
scatter(long, lat, c=z)
s1.set_aspect(1.)
s2 = subplot(122)
scatter(x, y, c=z)
s2.set_aspect(1.)