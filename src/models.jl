"""
Basic types and methods.

    DensityModel : parameterized density model
    AnisotropyModel : parameterized orbital anisotropy model
    JeansModel : collection of mass model, tracer density model, and anisotropy model
"""
module Models

using Slomo.Constants: G, rmax
using Slomo.Integrate: solve, integrate

struct NotImplemented <: Exception end

abstract type Model end

abstract type DensityModel <: Model end

abstract type AnisotropyModel <: Model end

struct JeansModel <: Model
    mass_model::Array{DensityModel, 1}
    density_model::DensityModel
    anisotropy_model::AnisotropyModel
end

JeansModel(mass_model::DensityModel,
           density_model::DensityModel,
           anisotropy_model::AnisotropyModel) = begin
               JeansModel([mass_model], density_model, anisotropy_model)
           end

"""
Create a new model from the specified parameters.
"""
function update(model::Model, parameters::Dict{Symbol, Float64})
    args = []
    for kw in propertynames(model)
        if kw in keys(parameters)
            push!(args, parameters[kw])
        else
            push!(args, getproperty(model, kw))
        end
    end
    return typeof(model)(args...)
end

update(model::Array{T} where T<:Model, parameters::Dict{Symbol, Float64}) = begin
    [update(m, parameters) for m in model]
end

update(model::JeansModel, parameters::Dict{Symbol, Float64}) = begin
    JeansModel(update(model.mass_model, parameters),
               update(model.density_model, parameters),
               update(model.anisotropy_model, parameters))
end

"""
Stupid hack to see if there is an analytic profile defined for a subtype
or if it needs to be numerically integrated.
"""
function has_analytic_profile(f, model::Model)
    method = methods(f, (typeof(model), Vararg)).ms[1]
    return occursin(string(typeof(model)), string(method))
end


#===========================
Functions for density models
===========================#

"""
Local volume density at r.
"""
function density(model::DensityModel, r)
    throw(NotImplemented("no defined density profile"))
end

"""
Local surface density at R.  If not defined for a subtype of DensityModel, then
calculate numerically as the Abel transform from the volume density profile.
"""
function density2d(model::DensityModel, R)
    integrand(r) = density(model, r) * r / √(r ^ 2 - R ^ 2)
    return 2 * integrate(integrand, R, rmax)
end

"""
Enclosed mass within r.  If not defined for a subtype of DensityModel, then
calculate numerically as the integral of the volume density.
"""
function mass(model::DensityModel, r)
    integrand(x) = 4pi * x ^ 2 * density(model, x)
    return integrate(integrand, eps() ^ (1/3), r)
end

mass(model::Array{T} where T<:DensityModel, r) = sum([mass(m, r) for m in model])

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
    throw(NotImplemented("no defined beta profile"))
end

"""
Integrating factor when solving for sigma_rad^2
"""
function g_jeans(model::AnisotropyModel, r)
    exp.(2.0 * integrate(x -> beta(model, x) / x, 1.0, r))
end

"""
Projection kernel for an anisotropy model.
"""
function K_jeans(model::AnisotropyModel, r, R)
    throw(NotImplemented("no defined Jeans Kernel"))
end

end
