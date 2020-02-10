"""
Collection of anisotropy models.

* [`IsotropicModel`](@ref): isotropic model
* [`ConstantBetaModel`](@ref): constant beta model
* [`RSBetaModel`](@ref): flexible beta model from Read & Steger
* [`beta`](@ref): orbital anisotropy profile
"""
module Anisotropy

import Slomo.Utils: B_inc, gamma
import Slomo.Models: Model
export IsotropicModel, ConstantBetaModel, RSBetaModel, beta

abstract type AnisotropyModel <: Model end

"""
    beta(model::AnisotropyModel, r)

Orbital anisotropy as a function of radius.  The orbital anisotropy is defined as 

```math
\\beta = 1 - \\sigma_\\theta^2 / \\sigma_r^2
```
"""
function beta(model::AnisotropyModel, r)
    throw(NotImplemented("no defined beta profile"))
end

"""
    g_jean(model::AnisotropyModel, r)

Integrating factor when solving for \$\\sigma_r^2\$
"""
function g_jeans(model::AnisotropyModel, r)
    exp.(2.0 * integrate(x -> beta(model, x) / x, 1.0, r))
end

"""
    K_jeans(model::AnisotropyModel, r, R)

Jeans projection kernel for an anisotropy model.
"""
function K_jeans(model::AnisotropyModel, r, R)
    throw(NotImplemented("no defined Jeans Kernel"))
end

"""
    IsotropicModel()

Parameter-less model representing an isotropic system (i.e. beta = 0).
"""
struct IsotropicModel <: AnisotropyModel end

beta(model::IsotropicModel, r) = zeros(size(r))

K_jeans(model::IsotropicModel, r, R) = begin
    @. sqrt(1.0 - 1.0 / (r / R) ^ 2)
end

g_jeans(model::IsotropicModel, r) = ones(size(r))


"""
    ConstantBetaModel(beta)

Model that assumes beta is constant with galactocentric radius.
"""
struct ConstantBetaModel <: AnisotropyModel
    beta::Float64
end
ConstantBetaModel() = ConstantBetaModel(0.0)

beta(model::ConstantBetaModel, r) = model.beta .* ones(size(r))

g_jeans(model::ConstantBetaModel, r) = r .^ (2 * model.beta)

K_jeans(model::ConstantBetaModel, r, R) = begin
    Γ = gamma
    β = model.beta
    # nudge half-integral beta
    if 2.0 * β - floor(2.0 * β) == 0.0
        β = β + rand([1, -1]) * eps()^0.5
    end
    u = r ./ R
    um2 = u .^ -2.0
    term1 = (1.5 - β) * √(π) * Γ(β - 0.5) / Γ(β)
    term2 = β * B_inc(β + 0.5, 0.5, um2)
    term3 = -B_inc(β - 0.5, 0.5, um2)
    return 0.5 * u .^ (2.0 * β - 1.0) .* (term1 .+ term2 .+ term3)
end


"""
    RSBetaModel(beta0, betaInf, rbeta, nbeta)

Flexible velocity anisotropy model from Read & Steger 2017 (eq. 9)

```math
\\beta(r) = \\beta_0 + (\\beta_\\infty - \\beta_0) / (1 + (r_\\beta / r))^{n_\\beta}
```

* `beta0`: inner asymptotic anisotropy
* `betaInf`: outer asymptotic anisotropy
* `rbeta`: transition radius
* `nbeta`: transition sharpness, higher is a faster transition

`beta0, betaInf, nbeta = 0, 1, 2` corresponds to the Osipkov-Merritt profile

`beta0, betaInf, nbeta = 0, 0.5, 1` corresponds to the Mamon-Lokas profile
"""
struct RSBetaModel <: AnisotropyModel
    beta0::Float64
    betaInf::Float64
    rbeta::Float64
    nbeta::Float64
end

RSBetaModel() = RSBetaModel(0.0, 0.5, 10.0, 1.0)

beta(model::RSBetaModel, r) = begin
    @. model.beta0 + (model.betaInf - model.beta0) / (1.0 + (model.rbeta / r) ^ model.nbeta)
end

g_jeans(model::RSBetaModel, r) = begin
    n = model.nbeta
    r_β = model.rbeta
    β_inf = model.betaInf
    Δβ = model.betaInf - model.beta0
    return @. r ^ (2.0 * β_inf) * ((r_β / r)^n + 1.0) ^ (2.0 * Δβ / n)
end

end
