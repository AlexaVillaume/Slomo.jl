# Slomo.jl Documentation

`Slomo.jl` is a Julia package for modeling the dynamics of spherically symmetric (e.g., slow-rotating) galaxies.

```@contents
```

## Getting started

Install Julia and `Slomo.jl`: [Installing Slomo](@ref).

```julia
using Slomo
mass_model = Halos.NFWModel()
tracer_model = SersicModel()
anisotropy_model = ConstantBetaModel()
jean_model = JeansModel(mass_model, tracer_model, anisotropy_model)
```

## Index

```@index
```
