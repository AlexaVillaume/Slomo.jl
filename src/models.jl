module Models

include("constants.jl")
include("integrate.jl")

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
Local volume density at r.
"""
function density(model::DensityModel, r)
    error("Must be called on subtype of DensityModel with defined volume density")
end

"""
Local surface density at R.  If not defined for a subtype of DensityModel, then
calculate numerically as the Abel transform from the volume density profile.
"""
function density2d(model::DensityModel, R)
    integrand_abel(r) = density(model, r) * r / √(r ^ 2 - Ri ^ 2)
    return 2 * [integrate(integrand_abel, Ri, Inf) for Ri in R]
end

"""
Enclosed mass within r.  If not defined for a subtype of DensityModel, then
calculate numerically as the integral of the volume density.
"""
function mass(model::DensityModel, r)
    integrand(r) = 4pi * r ^ 2 * density(model, r)
    return integrate(integrand, 0.0, r)
end

"""
Gravitational potential at r.  Equal to the square of the circular velocity.
If not defined for a subtype of DensityModel, then calculate from the enclosed
mass.
"""
function potential(model::DensityModel)
    return -G * mass(model, r) ./ r
end

#==============================
Functions for anisotropy models
==============================#

"""
Velocity anisotropy (beta) as a function of radius.  The velocity anisotropy is
defined as 1 - sigma_tan^2 / sigma_rad^2.
"""
function beta(model::AnisotropyModel, r)
    error("Must be called on subtype of AnisotropyModel with defined beta profile")
end

"""
Integrating factor when solving for sigma_rad^2
"""
function g_jeans(model::AnisotropyModel, r)
    exp.(2 * integrate(x -> beta(model, x) / x, 1.0, r))
end

"""
Projection kernel for an anisotropy model.
"""
function K_jeans(model::AnisotropyModel, r, R)
    error("Must be called on subtype of AnisotropyModel with defined Jeans kernel")
end

end
