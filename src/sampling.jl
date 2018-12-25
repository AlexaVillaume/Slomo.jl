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
using HDF5
using EllipsisNotation
using ProgressMeter: @showprogress


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
    res = pmap(y -> logp(y, kwargs...), [x[i, :] for i in 1:nwalkers])
    # deal with possible blobs
    if length(res[1]) > 1
        lp = [res[i][1] for i in 1:nwalkers]
        nblobs = length(res[1]) - 1
        # blobs is (nblobs, nwalker)-shaped array of blobs
        blobs = Array{Any, 2}(undef, nblobs, nwalkers)
        for i in 1:nblobs
            for k in 1:nwalkers
                blobs[i, k] = res[k][1 + i]
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
    new_blobs : additional returns from walkers, shape (nblobs, nwalkers, ..)
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
                    for b in 1:size(new_blobs)[1]
                        if length(size(blobs[b, k])) > 0
                            # blob b is an array
                            new_blobs[b, j][..] = blobs[b, k][..]
                        else
                            # blob b is a scalar
                            new_blobs[b, j] = blobs[b, k]
                        end
                    end
                end
            end
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

"""
Draw niter samples from logp and store them in fname.

If the file does not exist or overwrite = true, create a new hdf5 file.  Else
if append = true, append samples to the file.

If creating a new file, then the initial walker positions (x0) need to be set.

If the log probability function has additional returns (aka "blobs"), then store
these as specified by blob_scheme.  This should be an array of pairs  mapping names
(as strings) to tuples of (Type, size) where size is a tuple of integers that 
describe the size of the expected array in that blob.

E.g., if logp returns (-1234.5, [1.0, π, 2.7], 42), then the blob scheme would be
    
    ```
    blob_scheme = ["some_blob" => (Float64, (3,)), "more_blobs" => (Int, ())]
    ```
"""
function sample(logp, niter::Int, fname::AbstractString;
                overwrite = false, append = false,
                x0 = nothing, blob_scheme = nothing,
                gw_scale_a = 2.0, kwargs...)
    
    if isfile(fname) && !(overwrite || append)
        throw("$fname exists, do you want to overwrite or append?")
    elseif isfile(fname) && append
        # get the starting positions
        if x0 != nothing
            @warn("ignoring the user-passed starting position in favor of the most recent sample")
        end
        # get chain size and last walker positions, last lp values
        h5open(fname, "cw") do file
            nsamples, nwalkers, ndim = size(file["chain"])
            x0 = reshape(file["chain"][nsamples, 1:nwalkers, 1:ndim], (nwalkers, ndim))
            lp0 = reshape(file["lp"][nsamples, 1:nwalkers], (nwalkers,))
            if "blobs" in names(file)
                nblobs = length(names(file["blobs"]))
                blobs0 = Array{Any, 2}(undef, nblobs, nwalkers)
                for (j, (name, (blob_type, blob_size))) in enumerate(blob_scheme)
                    for k in 1:nwalkers
                        idx = map(x -> isa(x, Int) ? x : (1:x.indices.stop),
                                  to_indices(file["blobs/$name"], (nsamples, k, ..)))
                        b = reshape(file["blobs/$name"][idx...], blob_size)
                        if blob_size == ()
                            blobs0[j, k] = b[1]
                        else
                            blobs0[j, k] = b
                        end
                    end
                end
            else
                blobs0 = nothing
            end
        end
    else
        # either fname doesn't exist or we're overwriting it
        @assert x0 != nothing "need to provide starting walker positions"
        x0 = convert_parameter_array(x0)
        lp0, blobs0 = compute_lp(logp, x0; kwargs...)
        @assert all(isfinite.(lp0)) "initial positions returned non-finite numbers"
        if blobs0 != nothing
            @assert size(blobs0)[1] == length(blob_scheme) "blob scheme doesn't match output"
        end
        nwalkers, ndim = size(x0)
        nsamples = 0
        h5open(fname, "w") do file
            dset = d_create(file, "chain", Float64,
                            ((1, nwalkers, ndim), (-1, nwalkers, ndim)),
                            "chunk", (1, nwalkers, ndim))
            dset = d_create(file, "lp", Float64,
                            ((1, nwalkers), (-1, nwalkers)),
                            "chunk", (1, nwalkers))
            dset = d_create(file, "accepted", Int,
                            ((1, nwalkers), (-1, nwalkers)),
                            "chunk", (1, nwalkers))
            if blob_scheme != nothing
                for (name, (blob_type, blob_size)) in blob_scheme
                    dset = d_create(file, "blobs/$name", blob_type,
                                    ((1, nwalkers, blob_size...), (-1, nwalkers, blob_size...)),
                                    "chunk", (1, nwalkers, blob_size...))
                end
            end
        end
    end
    
    old_x = x0
    old_lp = lp0
    old_blobs = blobs0
    @showprogress 1 "Sampling..." for i in (nsamples + 1):(nsamples + niter)
        # draw sample
        new_x, new_lp, new_accepted, new_blobs = sample(logp, old_x;
                                                        lp0 = old_lp,
                                                        blobs0 = old_blobs,
                                                        gw_scale_a = gw_scale_a,
                                                        kwargs...)
        # save sample
        h5open(fname, "r+") do file
            
            # expand arrays
            set_dims!(file["chain"], (i, nwalkers, ndim))
            set_dims!(file["lp"], (i, nwalkers))
            set_dims!(file["accepted"], (i, nwalkers))
            
            # store chain, logp value, and accepted bits
            file["chain"][i, :, :] = new_x
            file["lp"][i, :] = new_lp
            file["accepted"][i, :] = new_accepted
            
            if blob_scheme != nothing
                for (j, (name, (blob_type, blob_size))) in enumerate(blob_scheme)
                    if blob_size == ()
                        # scalar blob
                        set_dims!(file["blobs/$name"], (i, nwalkers))
                        file["blobs/$name"][i, :] = new_blobs[j, :]
                    else
                        # array blob
                        set_dims!(file["blobs/$name"], (i, nwalkers, blob_size...))
                        for k in 1:nwalkers
                            # fudge to list of indices for arbitrary sized array
                            idx = map(x -> isa(x, Int) ? x : (1:x.indices.stop),
                                      to_indices(file["blobs/$name"], (i, k, ..)))
                            file["blobs/$name"][idx...] = new_blobs[j, k][..]
                        end
                    end
                end
            end
            old_x = new_x
            old_lp = new_lp
            old_blobs = new_blobs
        end
    end     
end
    
end    
