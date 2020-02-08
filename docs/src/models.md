# Models

## Overview

Models in `Slomo.jl` are specified as new types whose attributes track the parameters of the model.

For instance, we can make a new model representing an NFW halo as

```julia
halo = Slomo.Halos.NFWModel(25.3, 3.7e6)
```

As shown below, the `NFWModel` has two parameters, `rs` and `rhos`, representing the scale radius and scale density of the halo in units of kpc and solar masses per cubic kpc, respectively.

```@docs
Halos.NFWModel
```

We can then query the enclosed mass of this model at a particular radius with the `mass` function, e.g.,

```julia
mass(halo, 5)
```

## Module contents

```@autodocs
Modules = [Slomo.Models]
```

