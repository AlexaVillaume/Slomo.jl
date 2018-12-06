using GSL: sf_hyperg_2F1, sf_gamma

import Slomo.Models: AnisotropyModel, beta, K_jeans, g_jeans

"""
Redefined incomplete beta function in terms of hypergeometric function 2F1.
"""
B_inc(a, b, z) = begin
    F(z) = sf_hyperg_2F1(a, 1.0 - b, a + 1.0, z)
    return @. a ^ -1.0  * z ^ a * F(z)
end


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
ConstantBetaModel() = ConstantBetaModel(0.01)

beta(model::ConstantBetaModel, r) = model.beta .* ones(size(r))

g_jeans(model::ConstantBetaModel, r) = r .^ (2 * model.beta)

K_jeans(model::ConstantBetaModel, r, R) = begin
    F = sf_hyperg_2F1
    Γ = sf_gamma
    β = model.beta
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
"""
struct RSBetaModel <: AnisotropyModel
    beta0::Float64
    betaInf::Float64
    rbeta::Float64
    nbeta::Float64
end

beta(model::RSBetaModel, r) = begin
    @. model.beta0 + (model.betaInf - model.beta0) / (1.0 + (model.rbeta / r) ^ model.nbeta)
end

g_jeans(model::RSBetaModel, r) = begin
    n = model.nbeta
    r0 = model.rbeta
    betaInf = model.betaInf
    delta_beta = betaInf - model.beta0
    return @. r ^ (2.0 * betaInf) * ((r0 / r)^n + 1.0) ^ (2.0 * delta_beta / n)
end


