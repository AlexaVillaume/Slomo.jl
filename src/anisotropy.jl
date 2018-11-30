import Slomo.Models: AnisotropyModel, beta, K_jeans, g_jeans

struct IsotropicModel <: AnisotropyModel end

beta(model::IsotropicModel, r::Array{Float64, 1}) = zero(r)

K_jeans(model::IsotropicModel, r, R) = begin
    sqrt.(1.0 .- 1.0 ./ (r ./ R) .^ 2)
end

g_jeans(model::IsotropicModel, r) = one(r)

struct ConstantBetaModel <: AnisotropyModel
    beta::Float64
end

beta(model::ConstantBetaModel, r) = model.beta .* ones(r)

K_jeans(model::ConstantBetaModel, r, R) = begin
    error("need to define this")
end

