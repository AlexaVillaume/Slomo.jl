import Slomo.Models: AnisotropyModel, beta, K_jeans, g_jeans

struct IsotropicModel <: AnisotropyModel
end

beta(model::IsotropicModel, r::Array{Float64, 1}) = zero(r)

K_jeans(model::IsotropicModel, r::Array{Float64, 1}, R::Array{Float64, 1}) = begin
    sqrt.(1.0 .- 1.0 ./ (r ./ R) .^ 2)
end

g_jeans(model::IsotropicModel, r::Array{Float64, 1}) = 1.0

struct ConstantBetaModel <: AnisotropyModel
    beta::Float64
end

beta(model::ConstantBetaModel, r::Array{Float64, 1}) = model.beta .* ones(r)

K_jeans(model::ConstantBetaModel, r::Array{Float64, 1}, R::Array{Float64, 1}) = begin
    error("need to define this")
end

