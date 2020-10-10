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
    @show pids = [(i-1)*nchainspersounding+1:i*nchainspersounding for i = 1:length(s)]
    # set up dBzdt operators for sounding numbers in s, at pids, with fnames
end
