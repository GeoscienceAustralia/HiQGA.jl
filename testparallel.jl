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
              
