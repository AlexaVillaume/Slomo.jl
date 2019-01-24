using Slomo: IsotropicModel, ConstantBetaModel, Halos, SersicModel, sigma_los

using Test

mass_models = [
    Halos.NFWModel(),
    SersicModel(5.0, 1.0, 1e8),
    [Halos.NFWModel(), SersicModel(10.0, 4.0, 1e11)]
]

density_models = [
    SersicModel(5.0, 0.9),
    SersicModel(5.0, 4.0)
]

anisotropy_models = [
    IsotropicModel(),
    ConstantBetaModel(-1),
    ConstantBetaModel(0),
    ConstantBetaModel(0.5)
]

jeans_models = [JeansModel(mass_model, density_model, anisotropy_model)
                for mass_model in mass_models
                for density_model in density_models
                for anisotropy_model in anisotropy_models]

radii = [
    1,
    1.0,
    [1],
    [1.0],
    [2.0, 1.0, 3.0],
    [1e-2, 1e-1, 1e0, 1e1, 1e2],
    collect(range(1e-1, stop = 1e1, length = 100))
]

for model in jeans_models
    for R in radii
        s = zeros(size(R))
        try
            s = sigma_los(model, R)
        catch err
            if isa(err, DomainError)
                @warn(err, model=model, R=R)
                continue
            end
        end
        @test all(isfinite.(s))
        @test size(s) == size(R)
    end
end

