"""
Utilities for MCMC sampling.

Uses the Goodman & Weare affine-invariant ensemble sampling algorithm.
Heavily inspired by Dan Foreman-Mackey's emcee code.

Samples can be saved to an HDF5 with the structure:

    /chain    : (nsamples, nwalkers, ndim)-shaped array of samples
    /lp       : (nsamples, nwalkers)-shaped array of log probability values
    /accepted : (nsamples, nwalkers)-shaped array of accepted (1) or rejected (0)
    /blobs    : group containing (nsamples, nwalkers)-shaped arrays for each
                additional return from the logp function
"""
module Sampling

using Distributed: pmap

const allowed_blob_type = Union{T, Array{T}} where T<:Union{Real, AbstractString}

"""
Allow parameter arrays (nwalkers, ndim) to be passed in as an array of parameter
arrays of shape (ndim,).

If x is a single point in parameter space, initialize a gaussian ball of points
of shape (nwalkers, ndim) around that point.  Defaults to 128 walkers.
"""
function convert_parameter_array(x::Array;
                                 nwalkers = 128, initial_spread = 1e-3)
    if typeof(x) <: Array{T, 2} where T <: Real
        return convert(Array{Float64, 2}, x)
    elseif typeof(x) <: Array{T, 1} where T <: Real
        # intialize gaussian ball around the input point in parameter space
        ndim = length(x)
        x = [x + initial_spread * randn(ndim) for i in 1:nwalkers]
    end
    x = convert(Array{Float64, 2}, hcat(x...)')
    return x
end

"""
Compute the log probability for each walker, catching and collecting any other 
outputs of the log probability function.  Maps out the function calls to
available workers with Distributed.pmap.
"""
function compute_lp(logp, x::Array{Float64, 2}; kwargs...)
    nwalkers, ndim = size(x)
    res = pmap(y -> logp(y; kwargs...), [x[i, :] for i in 1:nwalkers])
    # deal with possible blobs
    if length(res[1]) > 1
        lp = [res[i][1] for i in 1:nwalkers]
        @assert(length(res[1]) <= 2,
                "logp should have at most two returns; the log probability and a dictionary of blobs")
        # blobs is (nwalker,)-shaped array of blob dictionaries
        blobs = Array{Dict{AbstractString, allowed_blob_type}, 1}(undef, nwalkers)
        for k in 1:nwalkers
            if res[k][2] != nothing
                blobs[k] = res[k][2]
            end
        end
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

Parameters

    logp : function: (params, kwargs...) -> float, blobs_dict
    x0 : walker positions, (nwalkers, ndim)-shaped array
    lp0 : log probability at x0, (nwalkers,)-shaped array,
        recalculate if not provided
    gw_scale_a : scale parameter in Goodman & Weare algorithm (defaults to 2)
    
Additional keyword arguments get passed to logp.

Returns

    new_x : sample of shape (nwalkers, ndim)
    new_lp : log probability at x1, shape (nwalkers,)
    new_accepted : whether or not proposal was accepted, shape (nwalkers,)
    new_blobs : array (nwalkers,) of dictionaries, with additional returns from logp
"""
function sample(logp, x0::Array{Float64, 2};
                lp0 = nothing, blobs0 = nothing,
                gw_scale_a = 2.0, kwargs...)
    a = gw_scale_a
    nwalkers, ndim = size(x0)
    if lp0 == nothing
        lp0, blobs0 = compute_lp(logp, x0; kwargs...)
    else
        @assert length(lp0) == nwalkers
    end
    new_x = copy(x0)
    new_lp = copy(lp0)
    new_accepted = falses(size(lp0))
    new_blobs = nothing
    if blobs0 != nothing
        new_blobs = copy(blobs0)
    end
    batch1 = Array{Int, 1}(1:(nwalkers / 2))
    batch2 = Array{Int, 1}((nwalkers / 2 + 1):nwalkers)
    divisions = [(batch1, batch2), (batch2, batch1)]
    for ensembles in divisions
        active, inactive = ensembles
        # compute proposed positons for active walkers
        u = rand(length(active))
        zs = @. ((a - 1.0) * u + 1.0) ^ 2 / a
        proposals = zeros((length(active), ndim))
        for i in 1:length(active)
            proposals[i, :] = zs[i] * x0[active[i], :] + (1.0 - zs[i]) * x0[rand(inactive), :]
        end
        # compute the log probability at the propsed walkers
        lp, blobs = compute_lp(logp, proposals)
        # get the walker indices of active walkers with finite logp values
        good_proposals = isfinite.(lp)
        # index into all walkers
        all_idx = active[good_proposals]
        # index into active walkers
        active_idx = (1:length(active))[good_proposals]
        for i in 1:length(all_idx)
            j = all_idx[i]
            k = active_idx[i]
            logratio = (ndim - 1) * log(zs[k]) + lp[k] - lp0[j]
            accepted = log(rand()) < logratio
            new_accepted[j] = accepted
            # save samples, logp values, and blobs for accepted samples
            if accepted
                new_x[j, :] = proposals[k, :]
                new_lp[j] = lp[k]
                if new_blobs != nothing
                    for name in keys(new_blobs[1])
                        new_blobs[j][name] = blobs[k][name]
                    end
                end
            end
        end
    end
    return new_x, new_lp, new_accepted, new_blobs
end
    
end    
