"""
Utilities for MCMC sampling.
"""
module Sampling

using Distributed: pmap
using ProgressMeter: @showprogress

"""
Allow parameter arrays (nwalkers, ndim) to be passed in as an array of parameter
arrays of shape (ndim,).
"""
function convert_parameter_array(x::Array)
    if typeof(x) <: Array{Float64, 2}
        return x
    end
    x = convert(Array{Float64, 2}, hcat(x...)')
    return x
end

"""
Compute the log probability for each walker, catching and collecting any other 
outputs of the log probability function.  Maps out the function calls to
available workers with Distributed.pmap.
"""
function compute_lp(logp, x::Array; kwargs...)
    x = convert_parameter_array(x)
    nwalkers, ndim = size(x)
    res = pmap(y -> logp(y, kwargs...), [x[i, :] for i in 1:nwalkers])
    # deal with possible blobs
    if length(res[1]) > 1
        lp = [res[i][1] for i in 1:nwalkers]
        blobs = [[res[i][2:end]...] for i in 1:nwalkers]
    else
        lp = res
        blobs = nothing
    end
    return lp, blobs
end


"""
Uses the ensemble sampling algorithm of Goodman & Weare to draw a new sample 
from the probability distribution, along with any blobs of data returned by the
probability function.

Heavily inspired by emcee (https://emcee.readthedocs.io/en/stable/) and the 
Julia implementation (https://github.com/madsjulia/AffineInvariantMCMC.jl).

Parameters

    logp : function: (params, kwargs...) -> float, blobs...
    x0 : walker positions, (nwalkers, ndim)-shaped array
    lp0 : log probability at x0, (nwalkers,)-shaped array,
        recalculate if not provided
    gw_scale_a : scale parameter in Goodman & Weare algorithm (defaults to 2)
    
Additional keyword arguments get passed to logp.

Returns

    new_x : sample of shape (nwalkers, ndim)
    new_lp : log probability at x1, shape (nwalkers,)
    new_accepted : whether or not proposal was accepted, shape (nwalkers,)
    blobs : additional returns from walkers, shape (nwalkers,)
"""
function sample(logp, x0::Array{Float64, 2};
                lp0 = nothing, gw_scale_a = 2.0, kwargs...)
    a = gw_scale_a
    nwalkers, ndim = size(x0)
    if lp0 == nothing
        lp0, blobs = compute_lp(logp, x0; kwargs...)
    else
        @assert length(lp0) == nwalkers
    end
    new_x = copy(x0)
    new_lp = copy(lp0)
    new_blobs = []
    new_accepted = zeros(nwalkers)
    batch1 = Array{Int, 1}(1:(nwalkers / 2))
    batch2 = Array{Int, 1}((nwalkers / 2 + 1):nwalkers)
    divisions = [(batch1, batch2), (batch2, batch1)]
    for ensembles in divisions
        active, inactive = ensembles
        u = rand(length(active))
        zs = @. ((a - 1.0) * u + 1.0) ^ 2 / a
        proposals = zeros((length(active), ndim))
        for i in 1:length(active)
            proposals[i, :] = zs[i] * x0[active[i], :] + (1.0 - zs[i]) * x0[rand(inactive), :]
        end
        lp, blobs = compute_lp(logp, proposals)
        logratio = @. (ndim - 1) * log(zs) + lp - lp0[active]
        accepted = log.(rand(length(active))) .< logratio
        new_x[active[accepted], :] = proposals[accepted, :]
        new_lp[active[accepted]] = lp[accepted]
        new_accepted[active] = accepted
        if blobs != nothing
            append!(new_blobs, blobs)
        end
    end
    return new_x, new_lp, new_accepted, new_blobs
end

function sample(logp, x0::Array, niter::Int;
                gw_scale_a = 2.0, kwargs...)
    x0 = convert_parameter_array(x0)
    nwalkers, ndim = size(x0)
    @assert nwalkers % 2 == 0 "Must have an even number of walkers!"
    @assert nwalkers > 2 * ndim "Need at least 2 walkers per dimension!"

    chain = NaN * ones((niter, nwalkers, ndim))
    lp = NaN * ones((niter, nwalkers))
    accepted = convert(Array{Int, 2}, zeros((niter, nwalkers)))
    blobs = []
    old_x = x0
    old_lp, old_blobs = compute_lp(logp, x0)
    @assert all(isfinite.(old_lp)) "One of the walkers returned an invalid probability!"
    
    @showprogress 1 "Sampling..." for i in 1:niter
        new_x, new_lp, new_accepted, new_blobs = sample(logp, old_x;
                                                        lp0 = old_lp,
                                                        gw_scale_a = gw_scale_a,
                                                        kwargs...)
        chain[i, :, :] = new_x
        lp[i, :] = new_lp
        accepted[i, :] = convert(Array{Int, 1}, new_accepted)
        append!(blobs, new_blobs)
        # set previous result
        old_x = new_x
        old_lp = new_lp
    end
    return chain, lp, accepted, blobs
end

end    
