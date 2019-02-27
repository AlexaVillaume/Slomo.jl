module Constants

using Cosmology: cosmology

const planck18 = cosmology(h = 0.6766,
                           Neff = 3.046,
                           OmegaM = 0.3111)                           

const planck15 = cosmology(h = 0.6774,
                           Neff = 3.046,
                           OmegaM = 0.3089)

const planck13 = cosmology(h = 0.6777,
                           Neff = 3.046,
                           OmegaM = 0.3071)

const wmap9 = cosmology(h = 0.6932,
                        Neff = 3.046,
                        OmegaM = 0.2865)

const default_cosmo = planck18

# Gravitational constant in Msun kpc (km / s)^2
const G = 4.30091727e-06

# critical density in cosmologically-invariant units of Msun h^2 / kpc3
const rho_crit_h2 = 2.774848e+02

# radians per arcsecond
const rad_per_arcsec = 4.84813681109536e-06

# maximum integration radius in kpc
const rmax = 1e2 


end
