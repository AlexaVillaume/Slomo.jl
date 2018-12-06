module Slomo

include("models.jl")
include("anisotropy.jl")
include("sersic.jl")
include("nfw.jl")
include("jeans.jl")

export NFWModel, SersicModel, ConstantBetaModel, IsotropicModel, RSBetaModel, Tracer,
    jeans, mass, density, density2d, potential, g_jeans, K_jeans, beta, update
    
end
