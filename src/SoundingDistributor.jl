module SoundingDistributor
using Distributed
export splittasks, getss, getpids

const nsafeworkers = 479

function splittasks(soundings::AbstractVector; nchainspersounding=nothing, ppn=nothing)
    nsoundings = length(soundings)
    ncores = nworkers()
    splittasks(;nsoundings, ncores, nchainspersounding, ppn)
end

function splittasks(;nsoundings=nothing, ncores=nothing, nchainspersounding=nothing, ppn=nothing)
    # split into sequential iterations of parallel soundings
    @assert !any(isnothing.([nsoundings, ncores, nchainspersounding, ppn]))
    @assert mod(ppn, nchainspersounding+1) == 0
    # crucial difference in letting global manager work or not
    subtractone = ncores > nsafeworkers ? 1 : 0
    nparallelsoundings = Int((ncores+1)/(nchainspersounding+1)) - subtractone
    nsequentialiters = ceil(Int, nsoundings/nparallelsoundings)
    @info "will require $nsequentialiters iterations of $nparallelsoundings soundings in one iteration"
    nsequentialiters, nparallelsoundings
end    

function getss(iter, nsequentialiters, nparallelsoundings, nsoundings)
    if iter<nsequentialiters
        ss = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
    else
        ss = (iter-1)*nparallelsoundings+1:nsoundings
    end
    ss
end

function getpids(i, nchainspersounding)
    if nworkers() <= nsafeworkers # make the global manager earn their pay
        pids = (i-1)*nchainspersounding+i:i*(nchainspersounding+1)
    else # prevent the lazy global manager getting overloaded
        pids = i*(nchainspersounding+1)+1:(i+1)*(nchainspersounding+1)
    end
    pids       
end

end