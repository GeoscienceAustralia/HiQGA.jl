using Distributed, LinearAlgebra
addprocs(2)
@everywhere using DistributedArrays, LinearAlgebra

function filldarray(a::AbstractArray)
    r = Array{Future, 1}(undef, nworkers())
    for (ipid, pid) in enumerate(workers())
        r[ipid] = @spawnat pid a
    end
    DArray(r)
end

function fillarray()
    r = [@spawnat pid myid()*ones(2) for pid in workers()]
    DArray(r)
end

function testclosure(d::DArray, a)
    function myclosure(a)
        localpart(d)[:] *= a
    end


    for(idx, pid) in enumerate(workers())
       remote_do(myclosure, pid, a)
    end
end

function fillarraywithnoarray()
    r = [@spawnat pid myid() for pid in workers()]
    DArray(r)
end

##
r = [@spawnat pid myid()*rand(100,100) for pid in workers()]
d = DArray(r, (200,100))

@everywhere fetchvalsvd(r::Future) = svd(fetch(r))

function dosvd(r::Array{Future, 1})
    s = Array{Any}(undef, nworkers())
    @sync for (i, rr) in enumerate(r)
        @async s[i] = remotecall_fetch(fetchvalsvd, rr.where, rr)
    end
    s
end

function dosvd(d::DArray)
    rr = Array{Future, 1}(undef, nworkers())
    @sync for (ipid, pid) in enumerate(procs(d))
        rr[ipid] = @spawnat pid svd(localpart(d))
    end
    fetch.(rr)
end

@everywhere fetchlocalsvd(d::DArray) = svd(localpart(d))

function dosvd2(d::DArray)
    s = Array{Any}(undef, nworkers())
    @sync for (ipid, pid) in enumerate(procs(d))
        @async s[ipid] = remotecall_fetch(fetchlocalsvd, pid, d)
    end
    s
end
