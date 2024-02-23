module RDP
using LinearAlgebra, Dates, ArchGDAL, Printf, 
    PyPlot, Images, FileIO, HiQGA
import GeoFormatTypes as GFT
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

function worldcoordinates(;gridr=nothing, gridz=nothing)
    # from https://support.esri.com/en-us/knowledge-base/faq-what-is-the-format-of-the-world-file-used-for-geore-000002860
    A = step(gridr)
    D = 0
    B = 0
    E = step(gridz)
    C = A/2
    F = gridz[1] + E/2
    [A, D, B, E, C, F]
end

makeplist(X,Y) = [[x,y] for (x,y) in zip(X,Y)]

reprojecttoGDA94(plist, fromepsg) = ArchGDAL.reproject(plist, GFT.EPSG(fromepsg), GFT.EPSG(epsg_GDA94) )

makexyfromlatlonglist(latlonglist) = [[l[i] for l in latlonglist] for i in 2:-1:1]

function writeworldXML(latlonglist; dst_dir= "", line="", suffix="", fnamebar="colorbar.jpg", gridz=nothing)
    img = load(joinpath(dst_dir, "jpeg", "$line"*suffix*".jpg"))
    h, w = size(img)
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
    ax.imshow(img; cmap, extent=[gridr[1], gridr[end], gridz[end], gridz[1]])
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

function doonecurtaintriad(line::Int; nlayers=0, pathname="", dr=0, dz=0, dst_dir="",
        donn=false, ϵfrac=0, src_epsg=0, barfigsize=(0.2, 1.2), isfirst=false, islast=false,
        dpi=400, cmap="turbo", fnamebar="colorbar.jpg", vmin=-Inf, vmax=Inf,
        writecolorbar=true, VE=20, shrink = 25_000)
    X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev = 
    transD_GP.readxyzrhoϕ(line, nlayers; pathname)
    plist = makeplist(X,Y)
    R = transD_GP.CommonToAll.cumulativelinedist(X,Y)
    pgood = rdpreduce(plist, ϵfrac*R[end])
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
        writeworldXML(latlonglist;line, gridz, suffix, dst_dir)
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
    nlayers=0, dr=0, dz=0, 
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
            barfigsize, dpi, cmap, fnamebar, writecolorbar=writecb, dst_dir, shrink, VE,
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
        transD_GP.writevtkfromxyzrho(ρlow, ρmid, ρhigh, X, Y, Z, zall, ln; dst_dir)
        transD_GP.writevtkphifromsummary(ϕmean, ϕsdev, X, Y, Z, ln; dst_dir)
        fn = joinpath(dst_dir, "Line_$(ln).vts")
        transD_GP.writevtkxmlforcurtain(fn; src_epsg=epsg_GDA94, dst_epsg=epsg_WGS84, suffix="", vmin, vmax)
    end
end

end # module