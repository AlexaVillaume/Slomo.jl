# Models

## Overview

Models in `Slomo.jl` are specified as new types whose attributes track the parameters of the model.

For instance, we can make a new model representing an NFW halo as

```julia
halo = Slomo.Halos.NFWModel(25.3, 3.7e6)
```

The `NFWModel` has two parameters, `rs` and `rhos`, representing the scale radius and scale density of the halo in units of kpc and solar masses per cubic kpc, respectively.  We can then query the enclosed mass of this model at a particular radius with the `mass` function, e.g.,

```julia
mass_at_5kpc = mass(halo, 5)
```

Models are immutable composite types (i.e. [`struct`](https://docs.julialang.org/en/v1/manual/types/index.html#Composite-Types-1)s in the usual Julia syntax).  For those familiar with Python, these work like [NamedTuples]().

While it isn't possible to change the value of a parameter in a `Model`, you can easily make a copy of the model with different parameters using the [`Slomo.Models.update`](@ref) function:

```julia
denser_halo = update(halo, rhos=4e6)
```

The `Models` module defines the basic functionality for working with `DensityModels` (a subtype of `Models`), including functions for evaluating the volume density, surface density, and enclosed mass profiles.  All subtypes of `DensityModel` should at least implement the `density()` method.  Most models implemented in Slomo already implement custom `mass` methods for computational convenience.

## Module contents

```@autodocs
Modules = [Slomo.Models]
```

