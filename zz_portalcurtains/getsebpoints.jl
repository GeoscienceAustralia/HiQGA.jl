using PyPlot
includet("/home/547/ar0754/.julia/dev/HiQGA/zz_portalcurtains/RDP.jl")
cd(@__DIR__)
dict = Dict("Xkey"=>"X", "Ykey"=>"Y", "elevkey"=>"ELEVATION", "segkey"=>"SegmentID")
X, Y, elev, seg = RDP.readegspoints("Menindee_Lakes_test_line_interp_MGA54_20240509.egs", dict)
## get xy corresponding to fine grid
using NearestNeighbors
R, gridr = transD_GP.getRandgridr(xrd, yrd, dr)
ncols = length(gridr)
xyfine, gridrfine = transD_GP.getallxyinr(xrd, yrd, 0.; rangelenth=ncols)
tree = KDTree(xyfine)
## get axp, axd first
# probabilistic
map(axp) do ax
    for a in ax[1:6]
        for (xx, yy, mahd) in zip(X, Y, elev)
            idxs, = nn(tree,[xx';yy'])
            a.plot(gridrfine[idxs], mahd, "--k", linewidth=2)
        end
    end
end
# deterministic
for a in axd[1:6]
    for (xx, yy, mahd) in zip(X, Y, elev)
        idxs, = nn(tree,[xx';yy'])
        a.plot(gridrfine[idxs], mahd, "--k", linewidth=2)
    end
end