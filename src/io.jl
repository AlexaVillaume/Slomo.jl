"""
HDF5 file I/O
"""
module IO

import Slomo.Sampling: sample, convert_parameter_array, compute_lp, allowed_blob_type

using HDF5
using EllipsisNotation
using ProgressMeter: @showprogress

export sample, init

"""
The state of one iteration of sampling.

    x : (nwalkers, ndim)-shaped array of walker positions
    lp : log probability, (nwalkers,)-shaped
    accepted : accepted (1) or not (0), (nwalkers,)-shaped
    blobs : array of blob dictionaries (one per walker)
"""
struct State
    x::Array{Float64, 2}
    lp::Array{Float64, 1}
    accepted::Array{Int, 1}
    blobs::Union{Array{Dict{AbstractString, allowed_blob_type}, 1}, Nothing}
end

State(x::Array{Float64, 2},
      lp::Array{Float64, 1},
      accepted::Array{Int, 1}) = State(x, lp, accepted, nothing)

State(x::Array{Float64, 2},
      lp::Array{Float64, 1};
      blobs::Union{Array{Dict{AbstractString, allowed_blob_type}, 1}, Nothing} = nothing) = begin
          accepted = ones(size(x)[1])
          State(x, lp, accepted, blobs)
      end

"""
Create a new hdf5 file.

    x0 : (nwalkers, ndim)-shaped array of initial walker positions
    blobs0 : (nwalkers,)-shaped array of blob dictionaries

"""
function init(x0::Array{Float64, 2}, hdf5_file::AbstractString;
              blobs0 = nothing, overwrite = false)
    
    mode = "cw"
    if isfile(hdf5_file) && !overwrite
        throw("$hdf5_file exists, do you want to overwrite it?")
    elseif overwrite
        mode = "w"
    end

    x0 = convert_parameter_array(x0)
    nwalkers, ndim = size(x0)
    
    # infer blob scheme
    if blobs0 != nothing
        blobs = blobs0[1]
        @assert(isa(blobs, Dict{S, allowed_blob_type} where S<:AbstractString),
                "blobs should be a dictionary mapping strings to numbers/string/arrays")
    else
        blobs = nothing
    end
    
    
    
    # create datasets for chain, lp, and accepted arrays
    h5open(hdf5_file, mode) do file
        dset = d_create(file, "chain", Float64,
                        ((1, nwalkers, ndim), (-1, nwalkers, ndim)),
                        "chunk", (1, nwalkers, ndim))
        dset = d_create(file, "lp", Float64,
                            ((1, nwalkers), (-1, nwalkers)),
                        "chunk", (1, nwalkers))
        dset = d_create(file, "accepted", Int,
                        ((1, nwalkers), (-1, nwalkers)),
                        "chunk", (1, nwalkers))
        if blobs != nothing
            for (name, blob) in blobs
                blob_size = size(blob)
                if blob_size == ()
                    # scalar blobs
                    blob_type = typeof(blob)
                else
                    # array blobs
                    blob_type = typeof(blob[1])
                end
                dset = d_create(file, "blobs/$name", blob_type,
                                ((1, nwalkers, blob_size...), (-1, nwalkers, blob_size...)),
                                "chunk", (1, nwalkers, blob_size...))
            end
        end
    end
    return hdf5_file
end

init(state::State, hdf5_file::AbstractString; overwrite = false) = init(state.x, hdf5_file;
                                                                        blobs0 = state.blobs,
                                                                        overwrite = overwrite)

"""
Append the i-th iteration to the hdf5 file.
"""
function append_state!(hdf5_file, state::State, i::Int)
    new_x = state.x
    new_lp = state.lp
    new_accepted = state.accepted
    new_blobs = state.blobs
    h5open(hdf5_file, "r+") do file
        nsamples, nwalkers, ndim = size(file["chain"])
        @assert i >= nsamples "don't want to overwrite data!"
        # expand arrays
        set_dims!(file["chain"], (i, nwalkers, ndim))
        set_dims!(file["lp"], (i, nwalkers))
        set_dims!(file["accepted"], (i, nwalkers))
        
        # store chain, logp value, and accepted bits
        file["chain"][i, :, :] = new_x
        file["lp"][i, :] = new_lp
        file["accepted"][i, :] = new_accepted

        if new_blobs != nothing
            blob_names = keys(new_blobs[1])
            for name in blob_names
                # expand blob array
                blob_size = size(new_blobs[1][name])
                set_dims!(file["blobs/$name"], (i, nwalkers, blob_size...))
                for k in 1:nwalkers
                    # fudge to list of indices for arbitrary sized hdf5 array
                    idx = map(x -> isa(x, Int) ? x : (1:x.indices.stop),
                              to_indices(file["blobs/$name"], (i, k, ..)))
                    file["blobs/$name"][idx...] = new_blobs[k][name]
                end
            end
        end
    end
end

"""
Retrieve the most recent iteration.
"""
function read_last_state(hdf5_file)
    # get last walker positions, last lp values, last accepted, and last blob states
    h5open(hdf5_file, "r") do file
        nsamples, nwalkers, ndim = size(file["chain"])
        x0 = reshape(file["chain"][nsamples, 1:nwalkers, 1:ndim], (nwalkers, ndim))
        lp0 = reshape(file["lp"][nsamples, 1:nwalkers], (nwalkers,))
        if "blobs" in names(file)
            blob_names = names(file["blobs"])
            blobs0 = Array{Dict{AbstractString, allowed_blob_type}, 1}(undef, nwalkers)
            for k in 1:nwalkers
                blobs = Dict{AbstractString, allowed_blob_type}
                for name in blob_names
                    # fudge to list of indices for arbitrary sized array
                    idx = map(x -> isa(x, Int) ? x : (1:x.indices.stop),
                              to_indices(file["blobs/$name"], (i, k, ..)))
                    blobs[name] = file["blobs/$name"][idx...]
                end
                blobs0[k] = blobs
            end
        else
            blobs0 = nothing
        end
    end
    return State(x0, lp0, accepted0, blobs0)
end

"""
Draw niter samples from logp and store them in hdf5_file.

    logp : log probability function: (params, kwargs...) -> float, blob_dict
"""
function sample(logp, niter::Int, hdf5_file::AbstractString;
                overwrite = false, append = false, x0 = nothing, gw_scale_a = 2.0,
                kwargs...)

    if isfile(hdf5_file) && !(overwrite || append)
        throw("$hdf5_file exists, do you want to overwrite or append?")
    elseif isfile(hdf5_file) && append
        # get the starting positions
        if x0 != nothing
            @warn("ignoring the user-passed starting position in favor of the most recent sample")
        end
        nsamples, nwalkers, ndim = size(f["chain"])
        state0 = read_last_state(hdf5_file)
        x0 = state0.x
        lp0 = state0.lp
        blobs0 = state0.blobs
    else
        # either hdf5_file doesn't exist or we're overwriting it, so create it
        @assert x0 != nothing "need to provide starting walker positions"
        x0 = convert_parameter_array(x0)
        nwalkers, ndim = size(x0)
        nsamples = 0
        lp0, blobs0 = compute_lp(logp, x0; kwargs...)
        @assert all(isfinite.(lp0)) "initial positions returned non-finite numbers"
        state0 = State(x0, lp0; blobs = blobs0)
        init(state0, hdf5_file; overwrite = overwrite)
    end

    # sample!
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
        state = State(new_x, new_lp, new_accepted, new_blobs)
        append_state!(hdf5_file, state, i)
        old_x = new_x
        old_lp = new_lp
        old_blobs = new_blobs
    end
end

end
