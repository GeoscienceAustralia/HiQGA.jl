using DelimitedFiles, PyPlot, NearestNeighbors, DataInterpolations
skipstart = 6
A = readdlm("VTEM-plus-7.3ms-pulse-darlingparoo.cfm"; skipstart)
##
f, ax = plt.subplots(3, 1, sharex=true, figsize=(7,8))
ax[1].plot(A[:,1], A[:,2])
ax[2].plot(0.5*(A[1:end-1,1]+A[2:end,1]), diff(A[:,2]))
idx = falses(size(A,1))
idx[1] = true
slopethreshfrac = 0.1
valuethreshfrac = 0.01
nwidth = 1
downstart = -0.001
δ = diff(A, dims=1)[:,2]
δ = [δ[1]; δ]
lastpointidx = 1
for ipoint = nwidth+1:length(idx)-nwidth
    global lastpointidx
    if A[ipoint,1] < downstart
        δavleft  = δ[ipoint-1]
        δavright = δ[ipoint+1]
        if (δavleft*δavright<0) 
            lastpointidx = ipoint
            idx[ipoint]  = true
        end    
        if ((idx[ipoint] & idx[ipoint-1]) & (abs((A[ipoint,2] - A[ipoint-1,2])/A[ipoint, 2]) < valuethreshfrac))
            idx[ipoint-1] = false
        end    
    else
        δlast = δ[lastpointidx]
        if abs((δ[ipoint]-δlast)/δ[ipoint-1]) > slopethreshfrac
            lastpointidx = ipoint
            idx[ipoint]  = true
        end        
    end
end
# also make sure first zero from shutoff is hit
idx[findfirst(isapprox.(A[:,2], 0.0; atol=1e-6))] = true
ax[1].plot(A[idx,1], A[idx,2],".-")
interp = LinearInterpolation(A[idx,2], A[idx,1])
winterp = interp.(A[:,1]) # Gives the linear interpolation value at all times
ax[3].plot(A[:,1], (winterp-A[:,2])./A[:,2]*100)
ax[3].set_xlabel("time s")
ax[3].set_ylabel("% difference")
ax[2].set_ylabel("slope")
ax[1].set_ylabel("amplitude")
ax[1].grid()
ax[2].grid()
ax[3].grid()
npoints = sum(idx)
ax[1].set_title("Slope change frac = $slopethreshfrac, number of segments = $(npoints-1)")
f.tight_layout()
open("hedgehog.jl", "w") do f
    write(f, "ramp = [\n")
end
io = open("hedgehog.jl", "a")     
writedlm(io, A[idx,:])
close(io)
open("hedgehog.jl", "a") do f
    write(f, "]\n")
end
