"""
Basic model types and methods
"""
module Models

using Slomo.Constants: G
using Slomo.Integrate: solve, integrate

export Model, DensityModel
export update, density, density2d, mass

struct NotImplemented <: Exception end

abstract type Model end

abstract type DensityModel <: Model end

"""
    update(model)

Create a new model from the specified parameters.
"""
function update(model::Model; parameters...)
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

update(model::Array{T} where T<:Model; parameters...) = begin
    [update(m; parameters...) for m in model]
end

"""
    has_analytic_profile(f, model::Model)

Stupid hack to see if there is an analytic profile defined for a subtype
or if it needs to be numerically integrated.
"""
function has_analytic_profile(f, model::Model)
    method = methods(f, (typeof(model), Vararg)).ms[1]
    return occursin(string(typeof(model)), string(method))
end

"""
    density(model::DensityModel, r)

Local volume density at r.
"""
function density(model::DensityModel, r)
    throw(NotImplemented("no defined density profile"))
end

"""
    density2d(model::DensityModel, R)

Local surface density at R.  If not defined for a subtype of `DensityModel`, then calculate
numerically as the Abel transform from the volume density profile.
"""
function density2d(model::DensityModel, R)
    integrand(r) = density(model, r) * r / √(r ^ 2 - R ^ 2)
    return 2 * integrate(integrand, R, rmax)
end

"""
    mass(model::DensityModel, r)

Enclosed mass within r.  If not defined for a subtype of `DensityModel`, then calculate 
numerically as the integral of the volume density.
"""
function mass(model::DensityModel, r)
    integrand(x) = 4π * x ^ 2 * density(model, x)
    return integrate(integrand, eps() ^ (1/3), r)
end

mass(model::Array{T} where T<:DensityModel, r) = sum([mass(m, r) for m in model])

"""
    potential(model::DensityModel, r)

Gravitational potential at r.  Equal to the square of the circular velocity.  If not 
defined for a subtype of `DensityModel`, then calculate from the enclosed mass.
"""
function potential(model::DensityModel, r)
    return -G * mass(model, r) ./ r
end

end
