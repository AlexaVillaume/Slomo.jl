using Slomo: mass, density, Halos, Constants
using Roots: ConvergenceFailed
using Test

rtol = 1e-3
maxevals = 100

# NFW
for Mvir in exp10.([8.0, 10.0, 12.0, 14.0])
    for mdef in ["vir", "200c", "200m", "500c"]
        for z in [0.0, 0.5, 1.0]
            if mdef == "vir" || mdef == "200c"
                cvir = Halos.hmcr(Mvir; mdef = mdef, z = z)
            else
                cvir = 10.0
            end
            halo = Halos.NFW_from_virial(Mvir, cvir; mdef = mdef, z = z)
            rs = Halos.scale_radius(halo)
            Rvir = Halos.virial_radius(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals)            
            @test isapprox(Mvir, mass(halo, Rvir), rtol = rtol)
            @test isapprox(cvir, Rvir / rs, rtol = rtol)
            @test isapprox(Mvir, Halos.virial_mass(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals), rtol = rtol)
            @test isapprox(cvir, Halos.concentration(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals), rtol = rtol)
            @test isapprox(Mvir, Halos.Mvir_from_Rvir(Rvir; mdef = mdef, z = z), rtol = rtol)
            @test isapprox(Rvir, Halos.Rvir_from_Mvir(Mvir; mdef = mdef, z = z), rtol = rtol)
            for rs_fraction in [1, 1.0, [1], [1.0], [0.5, 1.0, 1.0]]
                r = rs .* rs_fraction
                @test size(mass(halo, r)) == size(r)
                @test size(density(halo, r)) == size(r)
            end
        end
    end
end

# GNFW
for gamma in [eps()^0.5, 0.5, 1.0, 1.5, 2.0 - eps()^0.5]
    for Mvir in exp10.([8.0, 10.0, 12.0, 14.0])
        for mdef in ["vir", "200c", "200m", "500c"]
            for z in [0.0, 0.5, 1.0]
                if mdef == "vir" || mdef == "200c"
                    cvir = Halos.hmcr(Mvir; mdef = mdef, z = z)
                else
                    cvir = 10.0
                end
                halo = Halos.GNFW_from_virial(Mvir, cvir, gamma; mdef = mdef, z = z)
                rs = Halos.scale_radius(halo)
                Rvir = Halos.virial_radius(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals)            
                @test isapprox(Mvir, mass(halo, Rvir), rtol = rtol)
                @test isapprox(cvir, Rvir / rs, rtol = rtol)
                @test isapprox(Mvir, Halos.virial_mass(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals), rtol = rtol)
                @test isapprox(cvir, Halos.concentration(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals), rtol = rtol)
                @test isapprox(Mvir, Halos.Mvir_from_Rvir(Rvir; mdef = mdef, z = z), rtol = rtol)
                @test isapprox(Rvir, Halos.Rvir_from_Mvir(Mvir; mdef = mdef, z = z), rtol = rtol)
                for rs_fraction in [1, 1.0, [1], [1.0], [0.5, 1.0, 1.0]]
                    r = rs .* rs_fraction
                    @test size(mass(halo, r)) == size(r)
                    @test size(density(halo, r)) == size(r)
                end
            end
        end
    end
end

# CoreNFW
for Mvir in exp10.([8.0, 10.0, 12.0, 14.0])
    for mdef in ["vir", "200c", "200m", "500c"]
        for z in [0.0, 0.5, 1.0]
            if mdef == "vir" || mdef == "200c"
                cvir = Halos.hmcr(Mvir; mdef = mdef, z = z)
            else
                cvir = 10.0
            end
            Re = 10.0
            t_sf = 10.0
            halo = Halos.CoreNFW_from_virial(Mvir, cvir, Re, t_sf; mdef = mdef, z = z)
            rs = Halos.scale_radius(halo)
            Rvir = Halos.virial_radius(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals)            
            @test isapprox(Mvir, mass(halo, Rvir), rtol = rtol)
            @test isapprox(cvir, Rvir / rs, rtol = rtol)
            @test isapprox(Mvir, Halos.virial_mass(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals), rtol = rtol)
            @test isapprox(cvir, Halos.concentration(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals), rtol = rtol)
            @test isapprox(Mvir, Halos.Mvir_from_Rvir(Rvir; mdef = mdef, z = z), rtol = rtol)
            @test isapprox(Rvir, Halos.Rvir_from_Mvir(Mvir; mdef = mdef, z = z), rtol = rtol)
            for rs_fraction in [1, 1.0, [1], [1.0], [0.5, 1.0, 1.0]]
                r = rs .* rs_fraction
                @test size(mass(halo, r)) == size(r)
                @test size(density(halo, r)) == size(r)
            end
        end
    end
end

# Einasto
# for Mvir in exp10.([8.0, 10.0, 12.0, 14.0])
#     for mdef in ["vir", "200c", "200m", "500c"]
#         for z in [0.0, 0.5, 1.0]
#             if mdef == "vir" || mdef == "200c"
#                 cvir = Halos.hmcr(Mvir; mdef = mdef, z = z)
#             else
#                 cvir = 10.0
#             end
#             alpha = 0.16
#             halo = Halos.Einasto_from_virial(Mvir, cvir, alpha; mdef = mdef, z = z)
#             rs = Halos.scale_radius(halo)
#             Rvir = cvir * rs
#             Mvir_calc = Mvir
#             cvir_calc = cvir
#             try
#                 Rvir = Halos.virial_radius(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals)
#             catch err
#                 @warn(err, attempt="Rvir", Mvir=Mvir, cvir=cvir, mdef=mdef, z=z)
#             end
#             try
#                 Mvir_calc = Halos.virial_mass(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals)
#             catch err
#                 @warn(err, attempt="Mvir", Mvir=Mvir, cvir=cvir, mdef=mdef, z=z)
#             end
#             try
#                 cvir_calc = Halos.concentration(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals)
#             catch err
#                 @warn(err, attempt="cvir", Mvir=Mvir, cvir=cvir, mdef=mdef, z=z)
#             end
#             @test isapprox(Mvir, mass(halo, Rvir), rtol = rtol)
#             @test isapprox(cvir, Rvir / rs, rtol = rtol)
#             @test isapprox(Mvir, Mvir_calc, rtol = rtol)
#             @test isapprox(cvir, cvir_calc, rtol = rtol)
#             @test isapprox(Mvir, Halos.Mvir_from_Rvir(Rvir; mdef = mdef, z = z), rtol = rtol)
#             @test isapprox(Rvir, Halos.Rvir_from_Mvir(Mvir; mdef = mdef, z = z), rtol = rtol)
#             for rs_fraction in [1, 1.0, [1], [1.0], [0.5, 1.0, 1.0]]
#                 r = rs .* rs_fraction
#                 @test size(mass(halo, r)) == size(r)
#                 @test size(density(halo, r)) == size(r)
#             end
#         end
#     end
# end


# SolNFW
# for Mvir in exp10.([8.0, 10.0, 12.0, 14.0])
#     for mdef in ["vir", "200c", "200m", "500c"]
#         for z in [0.0, 0.5, 1.0]
#             if mdef == "vir" || mdef == "200c"
#                 cvir = Halos.hmcr(Mvir; mdef = mdef, z = z)
#             else
#                 cvir = 10.0
#             end
#             rs = 0.0
#             Rvir = 0.0
#             m22 = 5.0
#             halo = nothing
#             try
#                 halo = Halos.SolNFW_from_virial(Mvir, cvir, m22; mdef = mdef, z = z)
#             catch err
#                 if isa(err, ConvergenceFailed)
#                     @warn(err, halo="SolNFW", Mvir=Mvir, cvir=cvir, mdef=mdef, z=z)
#                     continue
#                 end
#             end
#             rs = Halos.scale_radius(halo)
#             Rvir = Halos.virial_radius(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals)            
#             @test isapprox(Mvir, mass(halo, Rvir), rtol = rtol)
#             @test isapprox(cvir, Rvir / rs, rtol = rtol)
#             @test isapprox(Mvir, Halos.virial_mass(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals), rtol = rtol)
#             @test isapprox(cvir, Halos.concentration(halo; mdef = mdef, z = z, rtol = rtol, maxevals = maxevals), rtol = rtol)
#             @test isapprox(Mvir, Halos.Mvir_from_Rvir(Rvir; mdef = mdef, z = z), rtol = rtol)
#             @test isapprox(Rvir, Halos.Rvir_from_Mvir(Mvir; mdef = mdef, z = z), rtol = rtol)
#             for rs_fraction in [1, 1.0, [1], [1.0], [0.5, 1.0, 1.0]]
#                 r = rs .* rs_fraction
#                 @test size(mass(halo, r)) == size(r)
#                 @test size(density(halo, r)) == size(r)
#             end
#         end
#     end
# end


