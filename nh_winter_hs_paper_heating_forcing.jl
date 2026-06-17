# NHWinterHeldSuarezPaperHeating modifies the original Held-Suarez thermal
# relaxation by keeping the tropospheric HS equilibrium temperature below
# 100 hPa, while replacing the stratospheric equilibrium temperature with a
# perpetual NH-winter polar-vortex profile. The stratospheric profile follows
# the Kushner and Polvani (2004) setup mirrored from SH winter to NH winter:
#
#   T_eq_strat(phi, p) = (1 - W_NH(phi)) T_US(p) + W_NH(phi) T_PV(p; gamma)
#
# where W_NH transitions around 50N with a 10 degree width, and gamma controls
# the winter polar stratospheric lapse rate. gamma_km = 4.0 corresponds to the
# colder / stronger-vortex case used by Kushner and Polvani.
#
# On top of that NH-winter background, this forcing adds the tropospheric
# diabatic heating perturbation used by Lindgren et al.:
#
#   Q'(lambda, phi, p) =
#       q0 * sin(m * (lambda - lambda0))
#          * exp(-0.5 * ((phi - phi0) / sigma_phi)^2)
#          * sin(pi * (-log(p / p0)) / log(pt / p0))
#       for pt <= p <= p0, and 0 otherwise.
#
# Default heating settings:
#   q0 = 6 K/day
#   m = 1
#   phi0 = 45N
#   sigma_phi = 0.175 radians, about 10.03 degrees
#   p0 = 800 hPa
#   pt = 200 hPa
#
# The forcing is zonal-mean-zero to leading order because of the sine factor.
# It is intended as an idealized tropospheric stationary-wave source with flat
# lower boundary, not as direct stratospheric heating.

using SpeedyWeather

Base.@kwdef mutable struct NHWinterHeldSuarezPaperHeating{NF, VectorType, MatrixType} <: SpeedyWeather.AbstractForcing
    sigma_boundary::NF = 0.7
    relax_time_slow::Second = Day(40)
    relax_time_fast::Second = Day(4)

    Tmin::NF = 200
    Tmax::NF = 315
    delta_T_y::NF = 60
    delta_theta_z::NF = 10

    tropopause_pressure_hpa::NF = 100
    gamma_km::NF = 4
    vortex_edge_latitude::NF = 50
    vortex_width_latitude::NF = 10

    q0_kday::NF = 6
    zonal_wavenumber::Int = 1
    heating_longitude_phase_deg::NF = 0
    heating_latitude_center_deg::NF = 45
    heating_latitude_sigma_deg::NF = 10.026761414789407
    heating_pressure_bottom_hpa::NF = 800
    heating_pressure_top_hpa::NF = 200

    us_tropopause_temp::NF = 216.65

    log_sigma::VectorType
    temp_relax_freq::MatrixType
    temp_equil_a::VectorType
    temp_equil_b::VectorType
    vortex_weight::VectorType
    wave_heating_horizontal::VectorType
end

function NHWinterHeldSuarezPaperHeating(spectral_grid::SpectralGrid; kwargs...)
    log_sigma = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        zeros(spectral_grid.NF, spectral_grid.nlayers),
    )
    temp_relax_freq = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        zeros(spectral_grid.NF, spectral_grid.nlayers, spectral_grid.nlat),
    )
    temp_equil_a = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        zeros(spectral_grid.NF, spectral_grid.nlat),
    )
    temp_equil_b = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        zeros(spectral_grid.NF, spectral_grid.nlat),
    )
    vortex_weight = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        zeros(spectral_grid.NF, spectral_grid.nlat),
    )
    wave_heating_horizontal = SpeedyWeather.on_architecture(
        spectral_grid.architecture,
        zeros(spectral_grid.NF, spectral_grid.npoints),
    )

    return NHWinterHeldSuarezPaperHeating{
        spectral_grid.NF,
        typeof(log_sigma),
        typeof(temp_relax_freq),
    }(
        ;
        log_sigma = log_sigma,
        temp_relax_freq = temp_relax_freq,
        temp_equil_a = temp_equil_a,
        temp_equil_b = temp_equil_b,
        vortex_weight = vortex_weight,
        wave_heating_horizontal = wave_heating_horizontal,
        kwargs...,
    )
end

function us_standard_temperature(p_pa)
    p = max(p_pa, 1.0)
    p20 = 5474.89
    temp20 = 216.65
    gas_constant = 287.0
    gravity = 9.80665

    if p >= p20
        return temp20
    end

    height_m = 20000.0 - gas_constant * temp20 / gravity * log(p / p20)
    return temp20 + 0.001 * (height_m - 20000.0)
end

function SpeedyWeather.initialize!(
    forcing::NHWinterHeldSuarezPaperHeating,
    model::SpeedyWeather.AbstractModel,
)
    (; coslat, sinlat, latd, londs, latds) = model.geometry
    sigma = model.geometry.σ_levels_full
    p0 = model.atmosphere.reference_pressure
    radius = model.planet.radius

    relax_slow = 1 / Second(forcing.relax_time_slow).value
    relax_fast = 1 / Second(forcing.relax_time_fast).value

    forcing.log_sigma .= log.(sigma)

    forcing.temp_relax_freq .= relax_slow .+
        (relax_fast - relax_slow) *
        max.(0, (sigma .- forcing.sigma_boundary) ./ (1 - forcing.sigma_boundary)) .*
        (coslat') .^ 4
    forcing.temp_relax_freq .*= radius

    forcing.temp_equil_a .= forcing.Tmax .-
        forcing.delta_T_y .* sinlat .^ 2 .+
        forcing.delta_theta_z .* log(p0) .* coslat .^ 2
    forcing.temp_equil_b .= -forcing.delta_theta_z .* coslat .^ 2

    forcing.vortex_weight .= 0.5 .* (
        1 .+
        tanh.((latd .- forcing.vortex_edge_latitude) ./ forcing.vortex_width_latitude)
    )

    tropopause_pressure_pa = 100 * forcing.tropopause_pressure_hpa
    forcing.us_tropopause_temp = us_standard_temperature(tropopause_pressure_pa)

    for ij in eachindex(forcing.wave_heating_horizontal)
        lat_factor = exp(
            -0.5 *
            ((latds[ij] - forcing.heating_latitude_center_deg) / forcing.heating_latitude_sigma_deg)^2,
        )
        wave_factor = sind(
            forcing.zonal_wavenumber *
            (londs[ij] - forcing.heating_longitude_phase_deg),
        )
        forcing.wave_heating_horizontal[ij] = lat_factor * wave_factor
    end

    return nothing
end

function SpeedyWeather.forcing!(
    vars::SpeedyWeather.Variables,
    forcing::NHWinterHeldSuarezPaperHeating,
    lf::Integer,
    model::SpeedyWeather.AbstractModel,
)
    temperature = vars.grid.temperature
    log_surface_pressure = vars.grid.pressure
    temperature_tendency = vars.tendencies.grid.temperature
    whichring = model.geometry.whichring

    kappa = model.atmosphere.κ
    reference_pressure = model.atmosphere.reference_pressure
    gas_constant = model.atmosphere.R_dry
    gravity = model.planet.gravity
    tropopause_pressure_pa = 100 * forcing.tropopause_pressure_hpa
    gamma = forcing.gamma_km / 1000
    polar_vortex_exponent = gas_constant * gamma / gravity
    heating_pressure_bottom_pa = 100 * forcing.heating_pressure_bottom_hpa
    heating_pressure_top_pa = 100 * forcing.heating_pressure_top_hpa
    heating_log_denominator = log(heating_pressure_top_pa / heating_pressure_bottom_pa)
    q0_scaled = forcing.q0_kday / 86400 * model.planet.radius

    for k in axes(temperature, 2)
        log_sigma_k = forcing.log_sigma[k]

        for ij in axes(temperature, 1)
            j = whichring[ij]
            log_pressure = log_surface_pressure[ij] + log_sigma_k
            pressure = exp(log_pressure)

            hs_equilibrium_temperature = max(
                forcing.Tmin,
                (
                    forcing.temp_equil_a[j] +
                    forcing.temp_equil_b[j] * log_pressure
                ) * (pressure / reference_pressure)^kappa,
            )

            if pressure < tropopause_pressure_pa
                us_temperature = us_standard_temperature(pressure)
                polar_vortex_temperature = forcing.us_tropopause_temp *
                    (pressure / tropopause_pressure_pa)^polar_vortex_exponent
                weight = forcing.vortex_weight[j]
                equilibrium_temperature = (1 - weight) * us_temperature +
                    weight * polar_vortex_temperature
            else
                equilibrium_temperature = hs_equilibrium_temperature
            end

            relaxation = forcing.temp_relax_freq[k, j]
            temperature_tendency[ij, k] -=
                relaxation * (temperature[ij, k] - equilibrium_temperature)

            if heating_pressure_top_pa <= pressure <= heating_pressure_bottom_pa
                vertical_factor = sin(
                    pi *
                    (-log(pressure / heating_pressure_bottom_pa)) /
                    heating_log_denominator,
                )
                temperature_tendency[ij, k] +=
                    q0_scaled * forcing.wave_heating_horizontal[ij] * vertical_factor
            end
        end
    end

    return nothing
end
