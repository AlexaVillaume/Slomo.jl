module Slomo

include("constants.jl")
include("integrate.jl")
include("models.jl")
include("anisotropy.jl")
include("sersic.jl")
include("nfw.jl")
include("jeans.jl")

using Slomo.Constants
using Slomo.Integrate
using Slomo.Models

export NFWModel, SersicModel, ConstantBetaModel, IsotropicModel, RSBetaModel, JeansModel,
    sigma_los, mass, density, density2d, potential, g_jeans, K_jeans, beta, update
    
end
