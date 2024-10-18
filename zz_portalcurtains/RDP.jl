module RDP
using LinearAlgebra, Dates, ArchGDAL, Printf, DataFrames, CSV,
    PyPlot, Images, FileIO, HiQGA, Interpolations, NearestNeighbors, DelimitedFiles,
    SegyIO
import GeoFormatTypes as GFT
import GeoDataFrames as GDF
import HiQGA.transD_GP.LineRegression.getsmoothline, HiQGA.transD_GP.CommonToAll.colstovtk

const mpl = PyPlot.matplotlib
const tilesize = 512
const epsg_GDA94 = 4283
const epsg_WGS84 = 4326

function perpdist(a, b, p)
    # a and b are two points on the line in space
    # p is the point to compute distance of from line
    n̂  = normalize(b-a)
    dperp = p - a - dot(p-a, n̂)n̂
    norm(dperp)
end

function rdpreduce(pointlist, ϵ=1.0)
    dists = perpdist.(Ref(pointlist[1]), Ref(pointlist[end]), pointlist)
    d, idx = findmax(dists)
    if d > ϵ
        # split line at the farthest point, doing RDP on left and right segments
        left  = rdpreduce(pointlist[1:idx], ϵ)
        right = rdpreduce(pointlist[idx:end], ϵ)
        # concatenate results of RDP on segments
        keep = vcat(left[1:end-1], right)
    else
        # don't keep the most distant point from the line, only the end points of the line
        keep = [pointlist[1], pointlist[end]]
    end    
    keep
end

function scaledRDP(xin...;ϵ=1.0) # input x,y,z,etc. 
    scaled = map(xin) do x
        (x .- minimum(x))/(maximum(x) - minimum(x))
    end
    plistscaled = map(zip(scaled...)) do (x)
         [x...]
    end
    scaledpgood = reduce(vcat, RDP.rdpreduce(plistscaled, ϵ)')
    scaledout = map(zip(eachcol(scaledpgood), xin)) do (s, x) 
        s*(maximum(x) - minimum(x)) .+ minimum(x)
    end
end

function worldcoordinates(;gridr=nothing, gridz=nothing)
    # from https://support.esri.com/en-us/knowledge-base/faq-what-is-the-format-of-the-world-file-used-for-geore-000002860
    A = step(gridr)
    D = 0
    B = 0
    E = step(gridz)
    C = A/2
    F = gridz[1] + E/2
    @sprintf("%.3f\n%.3f\n%.3f\n%.3f\n%.3f\n%.3f", A, D, B, E, C, F)
end

function makeextent(lnum, ncols, nrows, gridr, gridz; suffix="")
    @sprintf("%i%s\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f", 
        lnum, suffix, 0, 0, ncols, -nrows, gridr[1], maximum(gridz), gridr[end], minimum(gridz))
end    

function XYZ_zmid_gridtoSEGY(σ, X, Y, Z; dr=nothing, zall=nothing, dz=nothing, fname="line", dst_dir="", suffix="",
            nanval=-6 #=-6 is in log10=#,)
    img, gridr, gridz, topofine, R = transD_GP.makegrid(σ, X, Y, Z; 
            dr, zall, dz)
    img[isnan.(img)] .= nanval       
    segypath = dst_dir
    xymid, = transD_GP.getallxyinr(X, Y, dr)
    xm, ym, = map(i->xymid[i,:], 1:2)
    rm = transD_GP.CommonToAll.cumulativelinedist(xm, ym)
    topom = (interpolate((gridr,), topofine, Gridded(Linear())))(rm)
    block = SeisBlock(Float32.(img))
    set_header!(block, :dt, dz*1000) # ms
    set_header!(block, :dtOrig, dz*1000) # ms
    set_header!(block, :nsOrig, size(img,1))
    set_header!(block, :SourceX, round.(Int, xm))
    set_header!(block, :SourceY, round.(Int, ym))
    set_header!(block, :GroupX, round.(Int, xm))
    set_header!(block, :GroupY, round.(Int, ym))
    set_header!(block, :CDPTrace, Array(0:length(xm)-1))
    set_header!(block, :Tracenumber, Array(0:length(xm)-1))
    set_header!(block, :Inline3D, 1)
    set_header!(block, :Crossline3D, Array(1:length(xm)))
    set_header!(block, :ElevationScalar, round.(Int, topom))
    set_header!(block, :DelayRecordingTime, -round(Int, maximum(gridz))) # shift to topo as start
    fname = joinpath(segypath, fname*"_"*suffix)
    segy_write(fname*".segy", block)
end    

function writesegyfromxyzrhodir(nlayers::Int; src_dir="", src_epsg=0, dst_dir="", dr=nothing, dz=nothing)
    @assert !isnothing(src_epsg)
    lines = transD_GP.getprobabilisticlinesfromdirectory(src_dir)
    isdir(dst_dir) || mkpath(dst_dir)
    ioproj = open(joinpath(dst_dir, "0000_projection.txt"), "w")
    write(ioproj, "EPSG: $src_epsg")
    close(ioproj)
    map(lines) do ln
        @info "doing line $ln"
        fname = "line_$ln"
        X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev = transD_GP.readxyzrhoϕ(ln, nlayers; pathname=src_dir)
        [XYZ_zmid_gridtoSEGY(-ρ, X, Y, Z; dr, zall, dz, dst_dir, fname, suffix=str) for (ρ, str) in zip([ρlow, ρmid, ρhigh],["high", "mid", "low"])]
    end
end

function writepathextentfiles(line, nrows, ncols, gridr, topo, gridz, X, Y; suffix="", dst_dir="", src_epsg=0)
    geompath = joinpath(dst_dir,"geometry")
    isdir(geompath) || mkpath(geompath)
    f = open(joinpath(geompath, "$line"*suffix*".path.txt"), "w")
    Δr = (gridr[end]-gridr[1])/ncols
    xyfine, gridrfine = transD_GP.getallxyinr(X, Y, Δr; rangelenth=ncols+1)
    x, y = map(i->xyfine[i,:], 1:2)
    topofine = (interpolate((gridr,), topo, Gridded(Linear())))(gridrfine)
    plist = makeplist(x, y)
    latlonglist = reprojecttoGDA94(plist, src_epsg)
    long, lat = makexyfromlatlonglist(latlonglist)
    for i = 1:ncols
        write(f, @sprintf("%i%s\t%i\t%.3f\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",line, suffix, i, gridrfine[i], i, x[i], y[i], long[i], lat[i], topofine[i]))
    end  
    close(f)
    f = open(joinpath(geompath, "$line"*suffix*".extent.txt"), "w")
    write(f, makeextent(line, ncols, nrows, gridr, gridz; suffix))
    close(f)
    jpegpath = joinpath(dst_dir,"jpeg")
    f = open(joinpath(jpegpath, "$line"*suffix*".jgw"), "w")
    write(f, worldcoordinates(;gridr=gridrfine, gridz))
    close(f)
end    

makeplist(X,Y) = [[x,y] for (x,y) in zip(X,Y)]

reprojecttoGDA94(plist, fromepsg) = ArchGDAL.reproject(plist, GFT.EPSG(fromepsg), GFT.EPSG(epsg_GDA94) )

makexyfromlatlonglist(latlonglist) = [[l[i] for l in latlonglist] for i in 2:-1:1]

function colstovtk(cols::Dict, fname::String, src_epsg::Int; decfactor=1, hasthick=true, islog10=false, prefix="")
    # outputs in long, lat
    Xc, Yc, Zc, σc, thickc, linec = map(x->get(cols, x, 0), ["X", "Y", "Z", "cond", "thick", "line"])
    X, Y, Z, σ, thick, lines = transD_GP.readcols([Xc, Yc, Zc, σc, thickc, linec], fname; decfactor)
    plist = makeplist(X, Y)
    latlonglist = reprojecttoGDA94(plist, src_epsg)
    long, lat = makexyfromlatlonglist(latlonglist)
    colstovtk(long, lat, Z, σ, thick, lines, fname; hasthick, islog10, prefix)
end

function getimgsize(line, dst_dir, suffix)
    img = load(joinpath(dst_dir, "jpeg", "$line"*suffix*".jpg"))
    h, w = size(img)
end

function writeworldXML(latlonglist; h=0, w=0, dst_dir= "", line="", suffix="", fnamebar="colorbar.jpg", gridz=nothing)
    ntilelevels = gettilelevels(w,h,tilesize)
    firstpart = 
    """
    <?xml version="1.0" encoding="UTF-8" standalone="yes" ?>
    <Layer version="1" layerType="CurtainImageLayer">
        <DisplayName>$(string(line)*suffix)</DisplayName>
        <Legend>$fnamebar</Legend>
        <Service serviceName="DelegatorTileService">
            <URL>./</URL>
        </Service>
        <Delegates>
            <Delegate>LocalRequester</Delegate>
            <Delegate>TransparentColorTransformer(255,255,255,0.2)</Delegate>
            <Delegate>ResizeTransformer($tilesize,$tilesize)</Delegate>
        </Delegates>
        <LastUpdate>$(now())</LastUpdate>
        <DatasetName>$(string(line)*suffix)</DatasetName>
        <DataCacheName>AWS_CACHE_NAME/$(string(line)*suffix)</DataCacheName>
        <ImageFormat>image/jpg</ImageFormat>
        <FormatSuffix>.jpg</FormatSuffix>
        <AvailableImageFormats>
            <ImageFormat>image/jpg</ImageFormat>
        </AvailableImageFormats>
        <NumLevels count="$ntilelevels" numEmpty="0" />
        <TileSize>
            <Dimension width="$tilesize" height="$tilesize" />
        </TileSize>
        <FullSize>
            <Dimension width="$w" height="$h" />
        </FullSize>
        <Path>
    """
    middlepart = writelatlonglist(latlonglist)
    lastpart = 
    """    
        </Path>
        <CurtainTop>$(gridz[1] + step(gridz)/2)</CurtainTop>
        <CurtainBottom>$(gridz[end] - step(gridz)/2)</CurtainBottom>
        <FollowTerrain>false</FollowTerrain>
        <Subsegments>2</Subsegments>
        <UseTransparentTextures>true</UseTransparentTextures>
        <ForceLevelZeroLoads>true</ForceLevelZeroLoads>
        <RetainLevelZeroTiles>image/dds</RetainLevelZeroTiles>
        <UseMipMaps>true</UseMipMaps>
        <DetailHint>0.5</DetailHint>
    </Layer>
    """
    xmlpath = joinpath(dst_dir,"xml_and_tiles")
    isdir(xmlpath) || mkpath(xmlpath)
    f = open(joinpath(xmlpath, "$line"*suffix*".xml"), "w")
    map([firstpart, middlepart, lastpart]) do part
        write(f, part)
    end
    close(f)
end

function writelatlonglist(latlonglist)
    s = ""
    for p in latlonglist
        if p != last(latlonglist) 
            s*=@sprintf("\t\t<LatLon units=\"degrees\" latitude=\"%.6f\" longitude=\"%.6f\" />\n", p[1], p[2])
        else
            s*=@sprintf("\t\t<LatLon units=\"degrees\" latitude=\"%.6f\" longitude=\"%.6f\" />", p[1], p[2])
        end
    end    
    s        
end

function writedosbatchfile(jpegsavename::String; isfirst=false, islast=false)
    lname = basename(jpegsavename)
    dst_dir = dirname(dirname(jpegsavename))
    fname = joinpath(dst_dir, "tileall.bat")
    if isfirst
        f = open(fname, "w")
        write(f, "@echo off\n")
    else
        f = open(fname, "a")    
    end
    write(f, "call ribbon.bat -tilesize $tilesize -noLayerDef -source jpeg\\$lname -output xml_and_tiles\\\n")
    islast && write(f, "pause\n")
    close(f)
    nothing
end

function pltpoint(ax, pts)
    pt = reduce(hcat, pts)'
    ax.plot(pt[:,1], pt[:,2], "*-")
end 

function writeimageandcolorbar(img::Array, gridr, gridz, line::Int; cmap="turbo", dpi=300, dst_dir="",
        vmin=-Inf, vmax=Inf, fnamebar="colorbar.jpg", barfigsize=(1,6), suffix="", isfirst=false, islast=false,
        writecolorbar=false, VE=20, shrink=25_000)
    if (isinf(vmin) || isinf(vmax))
        vmin, vmax = extrema(filter(!isnan, img))
    end
    jpegpath = joinpath(dst_dir,"jpeg")
    if writecolorbar    
        fig, ax = plt.subplots(1, 2, gridspec_kw=Dict("width_ratios" => [0.3,1]), figsize=barfigsize)
        norm = mpl.colors.Normalize(vmin, vmax)
        cbar = fig.colorbar(mpl.cm.ScalarMappable(;norm, cmap), cax=ax[1], orientation="vertical")
        cbar.ax.tick_params(labelsize=4, pad=0)
        cbar.set_label("Log₁₀ S/m", fontsize=4, labelpad=0)
        ax[2].axis("off")
        fig.tight_layout()
        isdir(jpegpath) || mkpath(jpegpath)
        xmlpath = joinpath(dst_dir,"xml_and_tiles")
        savefig(joinpath(jpegpath, fnamebar); dpi)
        isdir(xmlpath) || mkpath(xmlpath)
        savefig(joinpath(xmlpath, fnamebar); dpi)
        close(fig)
    end
    figsize = gridr[end]/shrink, abs(diff([extrema(gridz)...])[1])*VE/shrink
    fig, ax = plt.subplots(;figsize)
    ax.imshow(img; cmap, extent=[gridr[1], gridr[end], gridz[end], gridz[1]], vmin, vmax)
    ax.set_aspect(VE)
    ax.set_aspect("auto")
    ax.axis("off")
    jpegsavename = joinpath(jpegpath, "$line"*suffix*".jpg")
    savefig(jpegsavename; dpi, bbox_inches = "tight", pad_inches = 0)
    writedosbatchfile(jpegsavename; isfirst, islast)
    close(fig)
    nothing
end

function gettilelevels(w,h,tilesize)
    xc = w/tilesize
    yc = h/tilesize
    levels = 0

    while (4 * xc * yc >= 1) 
        levels+=1;
        xc /= 2.
        yc /= 2.
    end
    levels
end    

function doonecurtaintriad(line::Int; nlayers=0, pathname="", dr=0, dz=0, dst_dir="", writegeom=false,
        donn=false, ϵfrac=0, src_epsg=0, barfigsize=(0.2, 1.2), isfirst=false, islast=false,
        dpi=400, cmap="turbo", fnamebar="colorbar.jpg", vmin=-Inf, vmax=Inf,
        writecolorbar=true, VE=20, shrink = 25_000)
    X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev = 
        transD_GP.readxyzrhoϕ(line, nlayers; pathname)
    X_, Y_, R_, = transD_GP.getXYlast(X, Y, dr)
    plist = makeplist(X_, Y_)
    pgood = rdpreduce(plist, ϵfrac*R_[end])
    latlonglist = RDP.reprojecttoGDA94(pgood, src_epsg)    
    map(zip([-ρhigh, -ρmid, -ρlow], ["_low", "_mid", "_high"],(1:3))) do (σ, suffix, i)    
        img, gridr, gridz, topofine, R = transD_GP.makegrid(σ, X, Y, Z; donn,
            dr, zall, dz)
        writecb = false
        if writecolorbar == true
            i == 1 && (writecb = true)
        end
        isl = false
        if islast == true
            i == 3 && (isl = true)
        end
        writeimageandcolorbar(img, gridr, gridz, line; cmap, fnamebar, suffix, dst_dir, isfirst=writecb, islast=isl,
        vmin, vmax, dpi, writecolorbar=writecb,
        barfigsize, VE, shrink)
        h, w = getimgsize(line, dst_dir, suffix)
        writeworldXML(latlonglist; h, w, line, gridz, suffix, dst_dir)
        writepathextentfiles(line, h, w, gridr, topofine, gridz, X_, Y_; suffix, dst_dir, src_epsg)
    end
end

function writebunchxmlfile(lines::Vector{Int}, dst_dir::String; prefix="", suffix="")
    xmlpath = joinpath(dst_dir,"xml_and_tiles")
    bunchfname = prefix*suffix
    f = open(joinpath(xmlpath, bunchfname*".xml"), "w")
    out = 
    """
    <DatasetList>
        <Dataset name="$bunchfname">
    """
    write(f, out)
    for line in lines
        write(f,"\t<Layer name=\"$(line)_$(suffix)\" url=\"$(line)_$(suffix).xml\"/>\n")
    end
    out = 
    """
        </Dataset>
    </DatasetList>
    """
    write(f, out)
    close(f)
    nothing
end    

function doallcurtaintriads(;src_dir="", dst_dir="curtains", prefix="",
    nlayers=0, dr=0, dz=0, writegeom=false,
    donn=false, ϵfrac=0, src_epsg=0, barfigsize=(0.2, 1.2), 
    dpi=400, cmap="turbo", fnamebar="colorbar.jpg", vmin=-Inf, vmax=Inf,
    VE=20, shrink = 25_000)
    
    isdir(dst_dir) || mkpath(dst_dir)
    lines = transD_GP.getprobabilisticlinesfromdirectory(src_dir)
    map(lines) do line
        @info "Doing line $line"
        writecb = line == first(lines) ? true : false
        isfirst = writecb
        islast = line == last(lines) ? true : false
        doonecurtaintriad(line; nlayers, pathname=src_dir, dr, dz, ϵfrac, src_epsg, isfirst, islast,
            barfigsize, dpi, cmap, fnamebar, writecolorbar=writecb, dst_dir, shrink, VE, writegeom,
            vmin, vmax)
    end
    for suffix in ("low", "mid", "high")
        writebunchxmlfile(lines, dst_dir; prefix, suffix)
    end    
    nothing
end

function writevtkfromxyzrhodir(nlayers::Int; src_dir="", dst_dir="", src_epsg=0, vmin=0, vmax=0)
    lines = transD_GP.getprobabilisticlinesfromdirectory(src_dir)
    isdir(dst_dir) || mkpath(dst_dir)
    map(lines) do ln
        X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev = transD_GP.readxyzrhoϕ(ln, nlayers; pathname=src_dir)
        plist = makeplist(X, Y)
        latlonglist = reprojecttoGDA94(plist, src_epsg)
        X, Y = makexyfromlatlonglist(latlonglist)
        transD_GP.CommonToAll.writevtkfromxyzrho(ρlow, ρmid, ρhigh, X, Y, Z, zall, ln; dst_dir)
        transD_GP.CommonToAll.writevtkphifromsummary(ϕmean, ϕsdev, X, Y, Z, ln; dst_dir)
        fn = joinpath(dst_dir, "Line_$(ln).vts")
        transD_GP.writevtkxmlforcurtain(fn; src_epsg=epsg_GDA94, dst_epsg=epsg_WGS84, suffix="", vmin, vmax)
    end
end

function writegiantfilefromxyzrhodir(nlayers::Int; src_dir="", src_epsg=0)
    lines = transD_GP.getprobabilisticlinesfromdirectory(src_dir)
    map(lines) do ln
        X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev = transD_GP.readxyzrhoϕ(ln, nlayers; pathname=src_dir)
        plist = makeplist(X, Y)
        latlonglist = reduce(hcat, reprojecttoGDA94(plist, src_epsg))'
        [latlonglist src_epsg*ones(size(X)) X Y Z transD_GP.zcentertoboundary(zall)'.*ones(size(X)) ρlow' ρmid' ρhigh' ϕmean ϕsdev]
    end
end

function writeshpfromxyzrhodir(nlayers::Int; prefix="", src_dir="", dst_dir="", src_epsg=0)
    lines = transD_GP.getprobabilisticlinesfromdirectory(src_dir)
    isdir(dst_dir) || mkpath(dst_dir)
    table = map(lines) do ln
        X, Y, _ = transD_GP.readxyzrhoϕ(ln, nlayers; pathname=src_dir)
        plist = makeplist(X, Y)
        longlat = reverse.(reprojecttoGDA94(plist, src_epsg))
        R = transD_GP.CommonToAll.cumulativelinedist(X,Y)
        (;geom=ArchGDAL.createlinestring(longlat), Line=ln, Length=R[end], soundings_inverted=length(R))
    end
    # writeesri
    path = joinpath(dst_dir, "shp")
    mkpath(path)
    fn = joinpath(path, prefix*".shp")
    GDF.write(fn, table; geom_columns=(:geom,), crs=GFT.EPSG(epsg_GDA94))
    # gpkg
    path = joinpath(dst_dir, "gpkg")
    mkpath(path)
    fn = joinpath(path, prefix*".gpkg")
    GDF.write(fn, table; geom_columns=(:geom,), crs=GFT.EPSG(epsg_GDA94))
end

function writeaseggdffromxyzrho(nlayers::Int; src_dir="", dst_dir="", 
         fname="", src_epsg=0, nunames=nothing, nuunits=nothing)
    isdir(dst_dir) || mkpath(dst_dir)
    sfmt = ["%15i", "%15.3f", "%15.3f", "%15.3f", "%15.3f", "%15.5f", "%15.5f", "%15.5f", "%15.5f", "%12.3f", "%12.3f"] 
    channel_names = [["Line", "X", "Y", "Z", "zcenter", "log10_cond_low", "log10_cond_mid", "log10_cond_high", "log10_cond_avg", 
                      "phid_mean", "phid_sdev"], 
                     ["", "m", "m", "m", "m", "Log10_Siemens_per_m", "Log10_Siemens_per_m", "Log10_Siemens_per_m", "Log10_Siemens_per_m",
                      "", ""],
                     ["Line", "X", "Y", "Z", "zcenter", "log10_cond_low", "log10_cond_mid", "log10_cond_high", "log10_cond_avg", 
                      "phid_mean", "phid_sdev"]
                    ]
    if !isnothing(nunames)
        nnu = length(nunames)
        sfmt = vcat(sfmt, fill("%15.4f", nnu*3)) # 3 times for lo, mid, hi as these are not vectors of same type
        channel_names[1] = vcat(channel_names[1], nunames.*"_low", nunames.*"_mid", nunames.*"_high")
        channel_names[3] = channel_names[1]
        channel_names[2] = vcat(channel_names[2], nuunits, nuunits, nuunits)
    end    
    outfile = joinpath(dst_dir, fname*"_EPSG_$src_epsg")
    lines = transD_GP.getprobabilisticlinesfromdirectory(src_dir)
    map(enumerate(lines)) do (iline, ln)
        @info "Doing Line $ln"
        if isnothing(nunames)
            X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev = transD_GP.readxyzrhoϕ(ln, nlayers; pathname=src_dir)
        else
            X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev, 
                                            nulow, numid, nuhigh = transD_GP.readxyzrhoϕnu(ln, nlayers; pathname=src_dir)
        end                                            
        for i in 1:length(X)
            mode = (iline == 1) & (i ==1) ? "w" : "a"
            ϕmeanwrite, ϕsdwrite = ϕmean[i], ϕsdev[i]
            if ϕmeanwrite > 1e4 # stop hogging the space in the column
                ϕmeanwrite, ϕsdwrite = 1e4, 1e4
            end
            vonerow = [ln, X[i], Y[i], Z[i], zall, -ρhigh[:,i], -ρmid[:,i], -ρlow[:,i], -ρavg[:,i], ϕmeanwrite, ϕsdwrite]
            if !isnothing(nunames)
                vonerow = vcat(vonerow, nulow[i,:], numid[i,:], nuhigh[i,:])
            end
            transD_GP.CommonToAll.writeasegdat(vonerow, sfmt, outfile, mode)
            if i == 1
                transD_GP.CommonToAll.writeasegdfnfromonerow(vonerow, channel_names, sfmt, outfile)
                transD_GP.dfn2hdr(outfile*".dfn")
            end    
        end    
    end
end

mutable struct XY
    x
    y
end 

function collectpoints(;npoints=1000)
    xy = ginput(npoints, timeout=0)
    x, y = [[pts[i] for pts in xy] for i =  1:2]
    XY(x, y)
end

function smoothline(xy::XY; λ²=0.01, finefactor=100, regtype=:R1, fname=nothing)
    xmin, xmax = extrema(xy.x)
    Δx = (xmax - xmin)/finefactor
    gridx = range(start=xmin, stop=xmax, step=Δx)
    gridy = snaptogrid(gridx, xy.x, xy.y)
    sd = 1 # identity weighting matrix for points
    ysmooth = getsmoothline(gridy, sd; δ²=λ², regtype)
    if !isnothing(fname)
        @assert !isfile(fname)
        io = open(fname,"w")
        write(io, "$(λ²)\n")
        write(io, "$finefactor\n")
        write(io, "$regtype\n")
        writedlm(io, [xy.x xy.y])
        close(io)
    end
    gridx, ysmooth
end    

function readpoints(fname::String)
    io = open(fname, "r")
    λ² = parse(Float64, readline(io))
    finefactor = parse(Float64, readline(io))
    regtype = Symbol(readline(io))
    xy = readdlm(io)
    xy = XY(xy[:,1], xy[:,2])
    smoothline(xy; λ², finefactor, regtype)
end

function readpoints(fnames::Vector{String})
    xy = map(fnames) do fn
        readpoints(fn)
    end
end

function snaptogrid(gridx, x, y)
    idx, _ = nn(KDTree(gridx'), x')
    gridy = NaN .+ zeros(size(gridx))
    gridy[idx] = y
    gridy
end

function readegspoints(fname::String, dict::Dict)
    musthave = ("Xkey", "Ykey", "elevkey", "segkey")
    map() do key
        @assert haskey(dict, key)
    end
    headers  = map(k->get(dict, k, 0), musthave)
    df = DataFrame(CSV.File(fname))
    out = map(headers) do h # X, Y, elev, seg
        df[:,h]
    end
    # split each segment
    segs = unique(out[end]) # the last one is segment ID
    X, Y, elev, seg = [[o[out[end] .== s] for s in segs] for o in out]
end

end # module