module Slomo

include("utils.jl")
include("constants.jl")
include("cosmology.jl")
include("integrate.jl")
include("models.jl")
include("halos/halos.jl")
include("sersic.jl")
include("anisotropy.jl")
include("jeans.jl")
include("sampling.jl")
include("io.jl")

using Slomo.Models
using Slomo.Halos
using Slomo.IO

export Halos, SersicModel, ConstantBetaModel, RSBetaModel, JeansModel
export sigma_los, mass, density, density2d, beta, update, sample

end
