#Steps to read and modify the SKyTEM headerfile 
using DataFrames
using CSV

#function to read the *.dfn file and its column number and write it as a *.hdr file 

function dfn2hdr(dfnfile::String,hdrfile::String)
    df = CSV.read(dfnfile, DataFrame; header=false);
    df[!, "Combined"] = string.(df[!, "Column3"], ":", df[!, "Column4"],":", df[!, "Column5"],":", df[!, "Column6"]);

    for i in 2:size(df,1)
        regex_literal = r"NAME="
        a1 = df[i,:].Column1[6:7]
        a2 = split(string(df[i,:].Combined), regex_literal)
        a2_r= first.(split.(a2[2], ":"))
        if occursin(";END DEFN", a2_r)
            a2_r = first.(split.(a2_r, ";"))
        end
        CSV.write(hdrfile,(Col_No = [a1], Col_Name = [a2_r]), append=true)
    end
end

dfn2hdr("GAB_10074_EM_Final.dfn","dfn2hdr.csv")