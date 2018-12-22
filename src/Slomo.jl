module Slomo

include("constants.jl")
include("cosmology.jl")
include("integrate.jl")
include("models.jl")
include("halos.jl")
include("sersic.jl")
include("anisotropy.jl")
include("jeans.jl")
include("sampling.jl")

using Slomo.Models
using Slomo.Halos
using Slomo.Sampling

export NFWModel, SersicModel, ConstantBetaModel, IsotropicModel, RSBetaModel, JeansModel,
    sigma_los, mass, density, density2d, potential, g_jeans, K_jeans, beta, update, sample

end
