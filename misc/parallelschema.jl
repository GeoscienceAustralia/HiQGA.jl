nsoundings = 24
ncores = 48
nchainspersounding = 4
@assert mod(ncores,nchainspersounding) == 0
nparallelsoundings = Int(ncores/nchainspersounding)
nsequentialiters = ceil(Int, nsoundings/nparallelsoundings)
@info "starting"
for iter = 1:nsequentialiters
    if iter<nsequentialiters
        @show s = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
    else
        @show s = (iter-1)*nparallelsoundings+1:nsoundings
    end
    @show pids = [(i-1)*nchainspersounding+2:i*nchainspersounding+1 for i = 1:length(s)]
    # set up dBzdt operators for sounding numbers in s, at pids, with fnames
end

## trying parallel, workers calling other workers
using Distributed
addprocs(4)
@info "workers are $(workers())"
##
@everywhere function sleeper(t)
    @info "sleeping at $(myid())"
    sleep(t)
end
@everywhere function trycall(t, pids)
    @sync for pid in pids
        @async remotecall_fetch(sleeper, pid, t)
    end
end
## now use them
ntimes = 2
tsleep = 2
@time for i = 1:ntimes
    @sync for (jpid, pids) in enumerate([[2,3],[4,5]])
        @async remotecall_fetch(trycall, jpid*2, tsleep, pids)
        println()
    end
end
