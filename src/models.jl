"""
Basic model types and methods
"""
module Models

using Slomo.Constants: G
using Slomo.Integrate: solve, integrate

export Model, DensityModel
export update, density, density2d, mass

struct NotImplemented <: Exception end

"""
    Model(parameters...)

A container for model parameters.  Model parameters can be aliased by passing in a 
dictionary of symbols mapping the aliased parameter name to the actual parameter name.

# Example

```julia
stars = SersicModel(10, 4, 1e6)
gcs = SersicModel(20, 3, 1, Dict(:Re_gcs => :Re))
```
"""
abstract type Model end

function Base.getproperty(model::Model, sym::Symbol)
    names = fieldnames(typeof(model))
    if sym in names
        return getfield(model, sym)
    elseif :aliases in names
        return getfield(model, getfield(model, :aliases)[sym])
    end
    getfield(model, sym)
end

function Base.propertynames(model::Model)
    Tuple(name for name in fieldnames(typeof(model)) if name != :aliases)
end

Base.show(io::IO, model::Model) = begin
    params = propertynames(model)
    values = Tuple(getproperty(model, p) for p in params)
    pstring = join(["$k=$v" for (k, v) in
                    zip(map(string, params), values)], ", ")
    print("$(typeof(model))($pstring)")
end

macro aliasable(expr)
    expr isa Expr && expr.head == :struct || error("Invalid usage of @aliasable")
    T = expr.args[2]
    if T isa Expr && T.head == :(<:)
        T = T.args[1]
    end
    params_ex = Expr(:parameters)
    call_args = Any[]
    Base._kwdef!(expr.args[3], params_ex.args, call_args)
    push!(expr.args[3].args, :(aliases::Dict{Symbol, Symbol}))
    constructor_expr = :(($(esc(T)))($(call_args...)) =
                         ($(esc(T)))($(call_args...), Dict{Symbol, Symbol}()))
    quote
        Base.@__doc__($(esc(expr)))
        $constructor_expr
    end
end

"""
    update(model)

Create a new model from the specified parameters.
"""
function update(model::Model; parameters...)
    args = []
    # aliases : aliased parameter name => true parameter name
    aliases = model.aliases
    inv_aliases = Dict(v => k for (k, v) in aliases)
    for name in propertynames(model)
        name = get(inv_aliases, name, name)
        if name in keys(parameters)
            push!(args, parameters[name])
        else
            push!(args, getproperty(model, name))
        end
    end
    push!(args, aliases)
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


#===============================================================================
DensityModel methods
===============================================================================#

abstract type DensityModel <: Model end

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
