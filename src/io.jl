"""
HDF5 file I/O
"""
module IO

using HDF5

struct HDF5Backend
    fname::AbstractString
end


"""
Create a new hdf5 file from a Julia source file.

The source file should define the following:

    nwalkers : Int, number of walkers to use
    logp : function (parameters, kwargs...) -> log probability
    x0 : (nwalker, ndim)-shaped array of initial walker positions

The logp function can also provide additional returns as arbitrary "blobs".
These should be passed as a dictionary mapping 
name::String => blob::Union{Array, Number}.

If the name of the output file is not specified, take from the source file,
replacing the *.jl extension with *.hdf5
"""
function init(source_file::AbstractString; hdf5_file = nothing, overwrite = false)
    
    @assert occursin(".jl", source_file) "source_file needs to be Julia code"
    if hdf5_file == nothing
        prefix = split(source_file, ".jl")[1]
        hdf5_file = prefix * ".hdf5"
    end
    if isfile(hdf5_file) && !overwrite
        throw("$hdf5_file exists, do you want to overwrite it?")
    end
    
    h5open(hdf5_file, "cw") do file
        source_lines = readlines(source_file)
        file["source"] = source_lines
    end

    # include the julia source code, super secure... hope for the best
    include(source_file)
    try
        x0
    catch err
        if isa(err, UndefVarError)
            throw("Need to define the initial walker positions as \"x0\"")
        else
            throw(err)
        end
    end
    try
        logp
    catch err
        if isa(err, UndefVarError)
            throw("Need to define the log probability function as \"logp\"")
        else
            throw(err)
        end
    end
    nwalkers, ndim = size(x0)
    
    # create datasets for chain, lp, and accepted arrays
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
    end    
    
    # infer blob scheme
    if length(res) > 1
        allowed_type = Union{Real,
                             AbstractString,
                             Array{T, N} where N where T<:Real,
                             Array{S, M} where M where S<:AbstractString}
        try
            res = logp(x0[1, :])
        catch err
            @warn "There was an error calling the logp function on one of the provided initial positions"
            throw(err)
        end
        @assert(isa(res[1], Real) && length(res) == 2,
                "logp should only return the log probability and a dictionary of blobs")
        blobs = res[2]
        
        @assert(isa(blobs, Dict{S, allowed_type} where S<:AbstractString),
                "blobs should be a dictionary mapping strings to numbers/string/arrays")
        h5open(hdf5_file, "cw") do file
            for (name, blob) in blobs
                blob_type = typeof(blob)
                blob_size = size(blob)
                dset = d_create(file, "blobs/$name", blob_type,
                                ((1, nwalkers, blob_size...), (-1, nwalkers, blob_size...)),
                                "chunk", (1, nwalkers, blob_size...))
            end
        end
    else
        blobs = nothing
    end    

    return HDF5Backend(hdf5_file)
end



end
