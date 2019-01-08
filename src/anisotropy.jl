import Slomo.Models: AnisotropyModel, beta, K_jeans, g_jeans
import Slomo.Utils: B_inc, gamma

#==============
Isotropic model
==============#

struct IsotropicModel <: AnisotropyModel end

beta(model::IsotropicModel, r) = zeros(size(r))

K_jeans(model::IsotropicModel, r, R) = begin
    @. sqrt(1.0 - 1.0 / (r / R) ^ 2)
end

g_jeans(model::IsotropicModel, r) = ones(size(r))

#==================
Constant beta model
==================#

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
        β = β + eps()
    end
    u = r ./ R
    um2 = u .^ -2.0
    term1 = (1.5 - β) * √(π) * Γ(β - 0.5) / Γ(β)
    term2 = β * B_inc(β + 0.5, 0.5, um2)
    term3 = -B_inc(β - 0.5, 0.5, um2)
    return 0.5 * u .^ (2.0 * β - 1.0) .* (term1 .+ term2 .+ term3)
end

#==================
Read & Steger model
==================#

"""
Flexible velocity anisotropy model from Read & Steger 2017 (eq. 9)

β(r) = β_0 + (β_inf - β_0) / (1 + (r_β / r)^n_β)
    
    beta0 : inner asymptotic anisotropy
    betaInf : outer asymptotic anisotropy
    rbeta : transition radius
    nbeta : transition sharpness, higher is a faster transition

beta0 = 0, betaInf = 1, nbeta = 2 is the Osipkov-Merritt profile

beta0 = 0, betaInf = 0.5, nbeta = 1 is the Mamon-Lokas profile
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


