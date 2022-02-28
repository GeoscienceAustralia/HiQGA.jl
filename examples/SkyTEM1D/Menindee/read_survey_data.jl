#Steps to read and modify the SKyTEM headerfile 
#module to read the hdr file and its column number 
using DataFrames
using CSV
using DelimitedFiles
using Statistics
using PyPlot

function read_survey_files(hdrfile::String,fname_dat::String,fname_specs_halt::String;
    #fname_dat= "",
    #fname_specs_halt="",
    frame_height = "",
    frame_dz = "",
    frame_dx = "",
    frame_dy = "",
    LM_Z = "",
    HM_Z = "",
    LM_Ã = "",
    HM_Ã = "",
    relerror = false,
    units=1e-12,
    figsize = (9,7),
    makesounding = false,
    dotillsounding = nothing,
    startfrom = 1,
    skipevery = 1,
    multnoise = 0.03,
    X = "",
    Y = "",
    Z = "",
    fid = "",
    linenum = "")
    df = CSV.read(hdrfile, DataFrame; header=false);
    ColNo,ColName = eachcol(df)
    # get column numbers for each required parameters 
    X = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("Easting", "easting"))), ColName)[1]),ColName)]))
    Y = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("Northing", "northing"))), ColName)[1]),ColName)]))
    Z = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("ground_elevation_lidar_aem_merged","GPS_Alt","gps_height"))),ColName)[1]),ColName)]))
    fid = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("Fiducial","fiducial"))),ColName)[1]),ColName)]))
    linenum = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("line", "Line", "Line number","line number","Line_Number","Line number"))),ColName)[1]),ColName)]))
    frame_height = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("frame_height","Frame_height"))),ColName)[1]),ColName)]))
    frame_dx = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("frame_dx","Frame_dx"))),ColName)[1]),ColName)]))
    frame_dy = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("frame_dy","Frame_dy"))),ColName)[1]),ColName)]))
    frame_dz = parse(Int64,(df[:,1][findfirst(==(filter(in(Set(("frame_dz","Frame_dz"))),ColName)[1]),ColName)]))
    LM_Z = parse.(Int64,split(df[:,1][findfirst(==(filter(in(Set(("LM_Z","Z_LM"))),ColName)[1]),ColName)],"-"))
    HM_Z = parse.(Int64,split(df[:,1][findfirst(==(filter(in(Set(("HM_Z","Z_HM"))),ColName)[1]),ColName)],"-"))

    # println("X = $X", "\n","Y = $Y", "\n" ,"Z = $Z", "\n" ,
    # "fid = $fid","\n" ,"linenum = $linenum","\n",
    # "frame_height = $frame_height","\n",
    # "frame_dx = $frame_dx","\n",
    # "frame_dy = $frame_dy","\n",
    # "frame_dz = $frame_dz","\n",
    # "LM_Z = $LM_Z","\n",
    # "HM_Z = $HM_Z","\n",
    # "relerror = false","\n",
    # "units = 1e-12")

    #call read_survey_file
    @assert frame_height > 0
    @assert frame_dz > 0
    @assert frame_dx > 0
    @assert frame_dy > 0
    @assert all(LM_Z .> 0)
    @assert all(HM_Z .> 0)
    if relerror
        @assert all(LM_Ã .> 0)
        @assert all(HM_Ã .> 0)
    end
    @assert X > 0
    @assert Y > 0
    @assert Z > 0
    @assert linenum > 0
    @assert fid > 0
    @info "reading $fname_dat"
    if dotillsounding!= nothing
        soundings = readdlm(fname_dat)[startfrom:skipevery:dotillsounding,:]
    else
        soundings = readdlm(fname_dat)[startfrom:skipevery:end,:]
    end
    easting = soundings[:,X]
    northing = soundings[:,Y]
    topo = soundings[:,Z]
    fiducial = soundings[:,fid]
    whichline = soundings[:,linenum]
    d_LM = soundings[:,LM_Z[1]:LM_Z[2]]
    d_HM = soundings[:,HM_Z[1]:HM_Z[2]]
    if relerror
        Ã_LM = soundings[:,LM_Ã[1]:LM_Ã[2]]
        Ã_HM = soundings[:,HM_Ã[1]:HM_Ã[2]]
    end
    zTx = soundings[:,frame_height]
    zRx = -(zTx + soundings[:,frame_dz])
    zTx = -zTx
    rRx = sqrt.(soundings[:,frame_dx].^2 + soundings[:,frame_dy].^2)

    @info "reading $fname_specs_halt"
    include(fname_specs_halt)
    @assert size(d_LM, 2) == length(LM_times)
    @assert size(d_HM, 2) == length(HM_times)
    if !relerror
        @assert size(d_LM, 2) == length(LM_noise)
        @assert size(d_HM, 2) == length(HM_noise)
        LM_noise[:] .*= units
        HM_noise[:] .*= units
    else
        Ã_LM[:] .*= units
        Ã_HM[:] .*= units
    end
    d_LM[:]     .*= units
    d_HM[:]     .*= units
    f = figure(figsize=figsize)
    ax = Array{Any, 1}(undef, 4)
    ax[1] = subplot(2,2,1)
    nsoundings = size(soundings, 1)
    plot_dLM = permutedims(d_LM)
    plot_dLM[plot_dLM .<0] .= NaN
    pcolormesh(1:nsoundings, LM_times, log10.(plot_dLM), shading="nearest")
    xlabel("sounding #")
    cbLM = colorbar()
    cbLM.set_label("log d_LM")
    ylabel("LM time s")
    axLM = ax[1].twiny()
    axLM.semilogy(LM_noise, LM_times)
    axLM.set_xlabel("high alt noise")
    ax[2] = subplot(2,2,3,sharex=ax[1], sharey=ax[1])
    plot_dHM = permutedims(d_HM)
    plot_dHM[plot_dHM .<0] .= NaN
    pcolormesh(1:nsoundings, HM_times, log10.(plot_dHM), shading="nearest")
    xlabel("sounding #")
    cbHM = colorbar()
    cbHM.set_label("log d_HM")
    ylabel("HM time s")
    ax[2].invert_yaxis()
    axHM = ax[2].twiny()
    axHM.semilogy(HM_noise, HM_times)
    axHM.set_xlabel("high alt noise")
    ax[3] = subplot(2,2,2, sharex=ax[1])
    plot(1:nsoundings, zRx, label="Rx")
    plot(1:nsoundings, zTx, label="Tx")
    legend()
    xlabel("sounding #")
    ylabel("height m")
    ax[3].invert_yaxis()
    ax[4] = subplot(2,2,4, sharex=ax[1])
    ax[4].plot(1:nsoundings, rRx)
    xlabel("sounding #")
    ylabel("rRx m")
    plt.tight_layout()
    if makesounding
        s_array = Array{SkyTEMsoundingData, 1}(undef, nsoundings)
        for is in 1:nsoundings
            l, f = Int(whichline[is]), fiducial[is]
            @info "read $is out of $nsoundings"
            dlow, dhigh = vec(d_LM[is,:]), vec(d_HM[is,:])
            if !relerror
                ÃLM = sqrt.((multnoise*dlow).^2 + LM_noise.^2)
                ÃHM = sqrt.((multnoise*dhigh).^2 + HM_noise.^2)
            else
                ÃLM, ÃHM = vec(Ã_LM[is,:]), vec(Ã_HM[is,:])
            end
            s_array[is] = SkyTEMsoundingData(rRx=rRx[is], zRxLM=zRx[is], zTxLM=zTx[is],
                zRxHM=zRx[is], zTxHM=zTx[is], rTx=rTx, lowpassfcs=lowpassfcs,
                LM_times=LM_times, LM_ramp=LM_ramp,
                HM_times=HM_times, HM_ramp=HM_ramp,
                LM_noise=ÃLM, HM_noise=ÃHM, LM_data=dlow, HM_data=dhigh,
                sounding_string="sounding_$(l)_$f",
                X=easting[is], Y=northing[is], Z=topo[is], fid=f,
                linenum=l)
        end
        return s_array
    end
end

read_survey_files("em_with_noise_calibration_lines.hdr",
"/g/data/z67/nfm547/transD_GP/examples/SkyTEM1D/Menindee/coincidewithtempest.dat",
"/g/data/z67/nfm547/transD_GP/examples/SkyTEM1D/Menindee/electronics_halt.jl")