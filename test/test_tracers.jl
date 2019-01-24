using Slomo: Models
using Slomo.Integrate: integrate

using Test

rtol = 1e-5

# sersic model
for model in [SersicModel(), SersicModel(1.0, 2.0), SersicModel(1.0, 2.0, 3.0)]
    Mtot = model.Mtot
    Re = model.Re
    integrand3d(x) = 4π * x^2 * density(model, x)
    integrand2d(x) = 2π * x * density2d(model, x)
    @test isapprox(Mtot, 2.0 * integrate(integrand2d, Re * rtol, Re), rtol = rtol)
    @test isapprox(Mtot, integrate(integrand2d, Re * rtol, Re / rtol), rtol = rtol)
    @test isapprox(Mtot, integrate(integrand3d, Re * rtol, Re / rtol), rtol = rtol)
    for R in [1, 1.0, [1], [1.0], [1.0, 2.0, 3.0]]
        @test size(density(model, R)) == size(R)
        @test size(density2d(model, R)) == size(R)
        @test size(mass(model, R)) == size(R)
        @test size(potential(model, R)) == size(R)
    end
end

# isotropic model
model = IsotropicModel()
for r in [1, 1.0, [1], [1.0], [1.0, 2.0, 3.0]]
    R = minimum(r) / 2.0
    @test beta(model, r) == zeros(size(r))
    @test size(beta(model, r)) == size(r)
    @test size(g_jeans(model, r)) == size(r)
    @test size(K_jeans(model, r, R)) == size(r)
end

# constant beta model
for model in [ConstantBetaModel(-30), ConstantBetaModel(-0.5), ConstantBetaModel(),
              ConstantBetaModel(0), ConstantBetaModel(0.5), ConstantBetaModel(0.96)]
    for r in [1, 1.0, [1], [1.0], [1.0, 2.0, 3.0]]
        R = minimum(r) / 2.0
        if isa(r, AbstractArray)
            β = beta(model, r)
            @test all(β .== β[1])
        end
        @test size(beta(model, r)) == size(r)
        @test size(g_jeans(model, r)) == size(r)
        @test size(K_jeans(model, r, R)) == size(r)
    end
end

# RS beta model
model = RSBetaModel()
for r in [1, 1.0, [1], [1.0], [1.0, 2.0, 3.0]]
    @test beta(model, r) == zeros(size(r))
    @test size(beta(model, r)) == size(r)
    @test size(g_jeans(model, r)) == size(r)
end

