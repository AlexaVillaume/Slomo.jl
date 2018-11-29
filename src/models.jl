module Models

abstract type Model end

abstract type DensityModel <: Model end

abstract type AnisotropyModel <: Model end

struct Tracer
    density_model::DensityModel
    anisotropy_model::AnisotropyModel
end

"""
Create a new model from the specified parameters.
"""

function update(model::Model, parameters::Dict{Symbol, Float64})
    args = []
    for kw in propertynames(model)
        try
            push!(args, parameters[kw])
        catch KeyError end
    end
    return typeof(model)(args...)
end

update(tracer::Tracer, parameters::Dict{Symbol, Float64}) = begin
    Tracer(update(tracer.density_model, parameters),
           update(tracer.anisotropy_model, parameters))
end

#===========================
Functions for density models
===========================#

"""
Local surface density at R.
"""
function surface_density(model::DensityModel, R::Array{Float64, 1})
    error("Must be called on subtype of DensityModel with defined surface density")
end

surface_density(model::DensityModel, R::Float64) = begin
    surface_density(model, collect([R]))
end

"""
Local volume density at r.
"""
function volume_density(model::DensityModel, r::Array{Float64, 1})
    error("Must be called on subtype of DensityModel with defined volume density")
end

volume_density(model::DensityModel, r::Float64) = begin
    volume_density(model, collect([r]))
end

"""
Enclosed mass within r.
"""
function mass(model::DensityModel, r::Array{Float64, 1})
    error("Must be called on subtype of DensityModel with defined mass")
end

mass(model::DensityModel, r::Float64) = begin
    mass(model, collect([r]))
end


#==============================
Functions for anisotropy models
==============================#

"""
Velocity anisotropy (beta) as a function of radius.  The velocity anisotropy is
defined as 1 - sigma_tan^2 / sigma_rad^2.
"""
function beta(model::AnisotropyModel, r::Array{Float64, 1})
    error("Must be called on subtype of AnisotropyModel with defined beta profile")
end

beta(model::AnisotropyModel, r::Float64) = beta(model, collect([r]))

"""
Projection kernel for an anisotropy model.
"""
function K_jeans(model::AnisotropyModel,
                 r::Array{Float64, 1},
                 R::Array{Float64, 1})
    error("Must be called on subtype of AnisotropyModel with defined beta profile")
end

K_jeans(model::AnisotropyModel, r::Float64, R::Float64) = begin
    return K_jeans(model, collect([r]), collect([R]))
end

"""
Integrating factor when solving for sigma_rad^2
"""
function g_jeans(model::AnisotropyModel, r::Array{Float64, 1})
    error("Must be called on subtype of AnisotropyModel with defined ")
end

g_jeans(model::AnisotropyModel, r::Float64) = g_jeans(model, collect([r]))

end
