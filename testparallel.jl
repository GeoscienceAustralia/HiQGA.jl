using Distributed
addprocs(2)
@everywhere using DistributedArrays

function filldarray(a::AbstractArray)
    r = Array{Future, 1}(undef, nworkers())
    for (ipid, pid) in enumerate(workers())
        r[ipid] = @spawnat pid a
    end
    DArray(r)
end    

