"""
    hmcr(Mvir; <keyword arguments>)

Compute the expected halo concentration from halo mass using the relation of 
Dutton & Maccio 2014 (equations 7, 10-13).

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
"""
function hmcr(Mvir; mdef = default_mdef, cosmo = default_cosmo, z = 0.0)
    mdef == lowercase(strip(mdef))
    if mdef == "200c"
        a = 0.520 + (0.905 - 0.520) * exp(-0.617 * z ^ 1.21)
        b = -0.101 + 0.026 * z
    elseif mdef == "vir"
        a = 0.537 + (1.025 - 0.537) * exp(-0.718 * z ^ 1.08)
        b = -0.097 + 0.024 * z
    else
        throw("only implemented for 200c and vir halo mass definitions")
    end
    # convert Mvir to h-scaled units
    return @. exp10(a + b * (log10(Mvir * cosmo.h) - 12.0))
end

"""
    hmcr_prior(Mvir, cvir; <keyword arguments>)
    
Compute the probability of drawing a halo with virial parameters (Mvir, cvir),
using the halo mass-concentration relation of Dutton & Maccio 2014.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `sigma_logc::Real = 0.16`: scatter in log_concentration at fixed Mvir (in dex)
"""
function hmcr_prior(Mvir, cvir;
                    mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                    sigma_logc = 0.16)
    logc_expected = log10(hmcr(Mvir; mdef = mdef, cosmo = cosmo, z = z))
    logc = log10(cvir)
    return log_gauss(logc, logc_expected, sigma_logc)
end

function hmcr_prior(halo::HaloModel;
                    mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                    sigma_logc = 0.16)
    Mvir = virial_mass(halo; mdef = mdef, z = z, cosmo = cosmo)
    logc_expected = log10(hmcr(Mvir; mdef = mdef, cosmo = cosmo, z =z))
    logc_expected = log10(cvir_expected)
    logc = log10(concentration(halo; mdef = mdef, z = z, cosmo = cosmo))
    return log_gauss(logc, logc_expected, sigma_logc)
end


"""
    shmr(Mvir; <keyword arguments>)

Compute the expected stellar mass from halo mass using the relation of 
Rodriguez-Puebla et al. 2017.

See equations 25-33 for the functional form and 49 - 55 for the parameters.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
"""
function shmr(Mvir; z = 0.0, cosmo = default_cosmo, mdef = default_mdef)

    if mdef != "vir"
        @warn("converting $mdef to vir halo mass definition", maxlog=1)
        @warn("assuming an NFW profile with standard concentration", maxlog=1)
        cvir = hmcr(Mvir; z = z, cosmo = cosmo, mdef = mdef)
        halo = NFW_from_virial(Mvir, cvir; z = z, cosmo = cosmo, mdef = mdef)
        Mvir = virial_mass(halo; z = z, cosmo = cosmo, mdef = "vir")
    end
    
    P(x, y, z) = y * z - x * z / (1.0 + z)
    Q(z) = exp(-4.0 / (1.0 + z)^2)
    
    ϵ0  =  -1.758
    ϵ1  =  +0.110
    ϵ2  =  -0.061
    ϵ3  =  -0.023
    M00 = +11.548
    M01 =  -1.297
    M02 =  -0.026
    α0  =  +1.975
    α1  =  +0.714
    α2  =  +0.042
    δ0  =  +3.390
    δ1  =  -0.472
    δ2  =  -0.931
    γ0  =  +0.498
    γ1  =  -0.157

    logϵ = ϵ0 + P(ϵ1, ϵ2, z) * Q(z) + P(ϵ3, 0.0, z)
    logM0 = M00 + P(M01, M02, z) * Q(z)
    α = α0 + P(α1, α2, z) * Q(z)
    δ = δ0 + P(δ1, δ2, z) * Q(z)
    γ = γ0 + P(γ1, 0.0, z) * Q(z)

    x = log10(Mvir) - logM0
    g(x) = δ * log10(1.0 + exp(x))^γ / (1.0 + exp(10.0^-x)) - log10(10^(-α * x) + 1.0)
    logMstar = logϵ + logM0 + g(x) - g(0)
    
    return exp10(logMstar)
end

function shmr(halo::HaloModel; mdef = default_mdef, z = 0.0, cosmo = default_cosmo)
    Mvir = virial_mass(halo; mdef = "vir", z = z, cosmo = cosmo)
    return shmr(Mvir; z = z, cosmo = cosmo, mdef = "vir")
end

"""
    shmr_prior(Mvir, Mstar; <keyword arguments>)

Compute the probability of drawing a halo mass and stellar mass pair (Mvir, Mstar).

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `sigma_logMstar::Real = 0.15` scatter in log_Mstar (in dex)
"""
function shmr_prior(Mvir, Mstar;
                    mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                    sigma_logMstar = 0.15)
    logMstar_expected = log10(shmr(Mvir; mdef = mdef, z = z, cosmo = cosmo))
    logMstar = log10(Mstar)
    return log_gauss(logMstar, logMstar_expected, sigma_logMstar)
end

function shmr_prior(halo::HaloModel, Mstar;
                    mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                    sigma_logMstar = 0.15)
    logMstar_expected = log10(shmr(halo; z = z, cosmo = cosmo))
    logMstar = log10(Mstar)
    return log_gauss(logMstar, logMstar_expected, sigma_logMstar)
end

"""
    peak_height(Mvir; <keyword arguments>)

Computes the dimensionless peak height,

```math
\\nu = \\delta_\\mathrm{crit}(z) / \\sigma(M, z)
```

using the approximation of Dutton & Maccio 2014 for a Planck 2013 cosmology.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
"""
function peak_height(Mvir; cosmo = default_cosmo, mdef="200c", z = 0.0)
    if mdef != "200c"
        throw("only implemented for 200c definition")
    end
    m = @. log10(M / 1e12 * cosmo.h)
    ν0 = @. exp10(-0.11 + 0.146 * m + 0.0138 * m^2 + 0.00123 * m^3)
    return @. ν0 * (0.033 + 0.79 * (1.0 + z) + 0.176 * exp(-1.356 * z))
end

"""
    alpha_peak(Mvir; <keyword arguments>)

Compute the Einasto shape parameter from peak height using the relation of Gao+2008.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
"""
function alpha_peak(Mvir; z = 0.0, cosmo = default_cosmo, mdef = "200c")
    if mdef != "200c"
        throw("only implemented for 200c definition")
    end
    ν = peak_height(Mvir, z; cosmo = cosmo, mdef = mdef)
    return @. 0.0095 * ν^2 + 0.155
end

"""
    alpha_peak_prior(Mvir, alpha; <keyword arguments>)

Compute the probability of drawing a halo with the given virial mass and Einasto shape
parameter using the alpha-peak height relation of Gao et al. 2008.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `sigma0::Real = 0.16`: scatter (in dex) of Mvir at redshift z = 0
- `sigmaz::Real = 0.03`: slope with redshift in scatter (in dex) of Mvir 
"""
function alpha_peak_prior(Mvir, alpha;
                          mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                          sigma0 = 0.16, sigmaz = 0.03)
    sigma_logalpha = sigma0 + sigmaz * z
    logalpha_expected = log10(alpha_peak(Mvir; z = z, cosmo = cosmo, mdef = mdef))
    logalpha = log10(alpha)
    return log_gauss(logalpha, logalpha_expected, sigma_logalpha)
end

"""
    abg_from_logshm(logshm)

αβγ scalings with log10(Mstar / Mhalo) from DiCintio+2014a.
"""
function abg_from_logshm(logshm)
    x = logshm
    if x < -4.1 || x > -1.3
        return 1, 3, 1
    end
    α = 2.94 - log10(exp10(-1.08 * (x + 2.33)) + exp10(2.29 * (x + 2.33)))
    β = 4.23 + 1.34 * x + 0.26 * x^2
    γ = -0.06 + log10(exp10(-0.68 * (x + 2.56)) + exp10(x + 2.56))
    return α, β, γ
end


