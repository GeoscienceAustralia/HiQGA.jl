module AEMwithNuisanceMcMCInversionTools
using ..AbstractOperator, ..CommonToAll
import ..AbstractOperator.makeoperator
import ..AbstractOperator.loopacrossAEMsoundings
import ..AbstractOperator.plotmodelfield!
import ..AbstractOperator.getndata
import ..Options, ..OptionsStat, ..OptionsNuisance
export makeAEMoperatorandnuisanceoptions, loopacrossAEMsoundings
import ..main # McMC function
using ..SoundingDistributor
using Distributed, Dates, Statistics, DelimitedFiles, PyPlot, Random

function make_tdgp_opt(sounding::Sounding;
                    rseed = nothing,
                    znall = [1],
                    fileprefix = "sounding",
                    nmin = 2,
                    nmax = 40,
                    K = GP.Mat32(),
                    demean = false,
                    sampledc = true,
                    sddc = 0.01,
                    sdpos = 0.05,
                    sdprop = 0.05,
                    fbounds = [-0.5 2.5],
                    λ = [2],
                    δ = 0.1,
                    pnorm = 2,
                    save_freq = 50,
                    nuisance_sdev   = [0.],
                    nuisance_bounds = [0. 0.],
                    updatenuisances = true,
                    C = nothing,
                    dispstatstoscreen = false,
                    restart = false
                    )
    sdev_dc = sddc*diff(fbounds, dims=2)[:]
    sdev_pos = [sdpos*abs(diff([extrema(znall)...])[1])]
    sdev_prop = sdprop*diff(fbounds, dims=2)[:]
    xall = permutedims(collect(znall))
    xbounds = permutedims([extrema(znall)...])

    history_mode = "w"
    restart && (history_mode = "a")

    updatenonstat = false
    needλ²fromlog = false
    if rseed != nothing
        Random.seed!(rseed)
    end

    opt = OptionsStat(fdataname = fileprefix*"_",
                            nmin = nmin,
                            nmax = nmax,
                            xbounds = xbounds,
                            fbounds = fbounds,
                            xall = xall,
                            λ = λ,
                            δ = δ,
                            demean = demean,
                            sampledc = sampledc,
                            sdev_dc = sdev_dc,
                            sdev_prop = sdev_prop,
                            sdev_pos = sdev_pos,
                            pnorm = pnorm,
                            quasimultid = false,
                            K = K,
                            save_freq = save_freq,
                            needλ²fromlog = needλ²fromlog,
                            updatenonstat = updatenonstat,
                            updatenuisances = updatenuisances,
                            dispstatstoscreen = dispstatstoscreen,
                            history_mode = history_mode
                            )

    bounds = makebounds(nuisance_bounds, sounding::Sounding)                   

    optn = OptionsNuisance(opt;
                    sdev = copy(nuisance_sdev),
                    bounds = bounds, C = C,
                    updatenuisances = updatenuisances)

                    opt, optn
                end
                
function makeAEMoperatorandnuisanceoptions(sounding::Sounding;
                                                    zfixed   = [-1e5],
                                                    ρfixed   = [1e12],
                                                    zstart = 0.0,
                                                    extendfrac = 1.06,
                                                    useML = false,
                                                    dz = 2.,
                                                    ρbg = 10,
                                                    nlayers = 40,
                                                    ntimesperdecade = 10,
                                                    nfreqsperdecade = 5,
                                                    showgeomplot = false,
                                                    plotfield = true,
                                                    nmin = 2,
                                                    nmax = 40,
                                                    K = GP.Mat32(),
                                                    demean = false,
                                                    sampledc = true,
                                                    sddc = 0.01,
                                                    sdpos = 0.05,
                                                    sdprop = 0.05,
                                                    fbounds = [-0.5 2.5],
                                                    λ = [2],
                                                    δ = 0.1,
                                                    save_freq = 50,
                                                    nuisance_sdev   = [0.],
                                                    nuisance_bounds = [0. 0.],
                                                    C = nothing,
                                                    updatenuisances = true,
                                                    dispstatstoscreen = false,
                                                    vectorsum = false
                                                    )
    
    aem, zall, znall, = makeoperator(sounding;
                                    zfixed, ρfixed, zstart, extendfrac,
                                    dz, ρbg, nlayers, ntimesperdecade,
                                    nfreqsperdecade, showgeomplot,
                                    plotfield, useML, vectorsum)
    
    opt, optn = make_tdgp_opt(sounding,
                                    znall = znall,
                                    fileprefix = soundings[idx].sounding_string,
                                    nmin = nmin,
                                    nmax = nmax,
                                    K = K,
                                    demean = demean,
                                    sampledc = sampledc,
                                    sddc = sddc,
                                    sdpos = sdpos,
                                    sdprop = sdprop,
                                    fbounds = fbounds,
                                    save_freq = save_freq,
                                    λ = λ,
                                    δ = δ,
                                    nuisance_bounds = nuisance_bounds,
                                    nuisance_sdev = nuisance_sdev,
                                    updatenuisances = updatenuisances,
                                    C = C,
                                    dispstatstoscreen = dispstatstoscreen
                                    )
    aem, opt, optn, zall
end

# Driver code for McMC inversion with no nuisances
function loopacrossAEMsoundings(soundings::Array{S, 1}, aem_in::Operator1D, opt_in::Options, optn_in::OptionsNuisance;
                            Tmax               = -1,
                            nsamples           = -1,
                            nchainsatone       =  1,
                            nchainspersounding = -1,
                            ppn                = -1) where S<:Sounding

    @assert ppn != -1
    @assert nchainspersounding != -1
    @assert nsamples != - 1
    @assert Tmax != -1

    nsoundings = length(soundings)
    nsequentialiters, nparallelsoundings = splittasks(soundings; nchainspersounding, ppn)

    for iter = 1:nsequentialiters
        ss = getss(iter, nsequentialiters, nparallelsoundings, nsoundings)
        @info "soundings in loop $iter of $nsequentialiters", ss
        @sync for (i, s) in Iterators.reverse(enumerate(ss))
            pids = getpids(i, nchainspersounding)
            @info "pids in sounding $s:", pids
            
            aem = makeoperator(aem_in, soundings[s])
            opt = deepcopy(opt_in)
            opt.fdataname = soundings[s].sounding_string*"_"
            optn = getoptnfromexisting(optn_in, opt, sounding)

            @async remotecall_wait(main, pids[1], opt, optn, aem, collect(pids[2:end]);
                                    Tmax, nsamples, nchainsatone)

        end # @sync
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end

end # module

