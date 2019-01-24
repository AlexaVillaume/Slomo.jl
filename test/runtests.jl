using Test

@time @testset "Cosmology" begin include("test_cosmo.jl") end
@time @testset "Tracer density and anisotropy models" begin include("test_tracers.jl") end
@time @testset "DM halo models" begin include("test_halos.jl") end
@time @testset "Jeans models" begin include("test_jeans.jl") end
# @time @testset "MCMC sampling" begin include("test_sampling.jl") end
# @time @testset "Utilities" begin include("test_utils.jl") end
