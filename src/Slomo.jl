module Slomo

include("models.jl")
include("anisotropy.jl")
include("sersic.jl")
include("nfw.jl")
include("jeans.jl")

export NFWModel, SersicModel, ConstantAnisotropyModel, IsotropicModel, jeans

end # module
