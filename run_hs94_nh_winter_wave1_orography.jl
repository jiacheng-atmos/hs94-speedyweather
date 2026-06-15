# This experiment starts from the original local HS94 baseline but adds the
# following modifications for an idealized SSW / stationary-wave experiment:
#
# 1. The model remains dry primitive-equation SpeedyWeather with Held-Suarez
#    style Newtonian temperature relaxation and linear lower-boundary drag.
# 2. The default HS94 equilibrium temperature is replaced by NHWinterHeldSuarez:
#    the troposphere keeps the HS profile, while p < 100 hPa uses a perpetual
#    NH-winter polar-vortex equilibrium temperature based on Kushner and
#    Polvani (2004), mirrored from their SH-winter setup to NH winter.
# 3. The stratospheric vortex strength is set by gamma_km = 4.0 K/km, the cold /
#    strong-vortex case in Kushner and Polvani's notation.
# 4. The vertical grid is L40 but no longer the SpeedyWeather default equally
#    spaced sigma grid. A custom nonuniform SigmaCoordinates grid puts the top
#    full level at sigma = 0.0005, which is about 0.5 hPa if ps = 1000 hPa.
# 5. A stationary wave-1 orography is added from 30N to 90N, with a sine
#    meridional envelope that is zero at 30N and 90N and peaks near 60N.
# 6. The time step is reduced from the default T31 value of 40 minutes to
#    20 minutes because the 0.5 hPa high-top custom sigma grid is less stable
#    than the original low-top HS94-style grid.
# 7. The output folder name records the key features: NH winter, gamma4,
#    wave-1 orography, T31L40, 0.5 hPa model top, and 20-day run length.
#
# This is not a literal reproduction of Kushner and Polvani (2004), because
# their controlled experiment used a flat lower boundary and varied gamma.
# Here gamma is fixed and the stationary wave forcing is varied through the
# imposed topography.

using SpeedyWeather

include("nh_winter_hs_forcing.jl")

Base.@kwdef mutable struct PrintDayProgress <: SpeedyWeather.AbstractCallback
    interval_days::Float64 = 1.0
    next_day::Float64 = 0.0
    total_days::Float64 = 0.0
end

function SpeedyWeather.initialize!(
    callback::PrintDayProgress,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    callback.next_day = callback.interval_days
    callback.total_days = Second(vars.prognostic.clock.period).value / 86400
    println("day 0.0 / ", callback.total_days)
    return nothing
end

function SpeedyWeather.callback!(
    callback::PrintDayProgress,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    elapsed_days = Millisecond(vars.prognostic.clock.time - vars.prognostic.clock.start).value / 86400000

    while elapsed_days + 1.0e-9 >= callback.next_day
        println("day ", callback.next_day, " / ", callback.total_days)
        callback.next_day += callback.interval_days
    end

    return nothing
end

SpeedyWeather.finalize!(::PrintDayProgress, args...) = nothing

truncation = 31
nlayers = 40
run_period = Day(20)
time_step_at_T31 = Minute(20)

sigma_half = [
    0.0,
    0.001,
    0.002,
    0.003,
    0.005,
    0.007,
    0.010,
    0.015,
    0.020,
    0.030,
    0.040,
    0.055,
    0.070,
    0.090,
    0.110,
    0.135,
    0.160,
    0.190,
    0.220,
    0.255,
    0.290,
    0.330,
    0.370,
    0.410,
    0.455,
    0.500,
    0.545,
    0.590,
    0.635,
    0.680,
    0.725,
    0.765,
    0.805,
    0.845,
    0.880,
    0.910,
    0.935,
    0.955,
    0.975,
    0.990,
    1.0,
]

gamma_km = 4.0
tropopause_pressure_hpa = 100.0
vortex_edge_latitude = 50.0
vortex_width_latitude = 10.0

orography_amplitude_m = 1000.0
zonal_wavenumber = 1
lat_south = 30.0
lat_north = 90.0

output_dir_id = "hs94_nh_winter_gamma4_wave1_orography_T31L40_top0p5hPa_dt20min_20d"

function wave1_orography(longitude, latitude)
    if lat_south <= latitude <= lat_north
        meridional_envelope = sind(180 * (latitude - lat_south) / (lat_north - lat_south))
        return orography_amplitude_m * cosd(zonal_wavenumber * longitude) * meridional_envelope
    end

    return 0.0
end

spectral_grid = SpectralGrid(
    trunc = truncation,
    nlayers = nlayers,
)
vertical_coordinates = SigmaCoordinates(sigma_half)
geometry = Geometry(
    spectral_grid;
    vertical_coordinates = vertical_coordinates,
)

forcing = NHWinterHeldSuarez(
    spectral_grid;
    gamma_km = gamma_km,
    tropopause_pressure_hpa = tropopause_pressure_hpa,
    vortex_edge_latitude = vortex_edge_latitude,
    vortex_width_latitude = vortex_width_latitude,
)
time_stepping = Leapfrog(
    spectral_grid;
    Δt_at_T31 = time_step_at_T31,
)
drag = LinearDrag(spectral_grid)
orography = ManualOrography(spectral_grid)
callbacks = CallbackDict(:progress => PrintDayProgress(interval_days = 1.0))
output = NetCDFOutput(
    spectral_grid,
    PrimitiveDry;
    id = output_dir_id,
    interval = Day(1),
)

model = PrimitiveDryModel(
    spectral_grid;
    geometry = geometry,
    forcing = forcing,
    drag = drag,
    orography = orography,
    time_stepping = time_stepping,
    callbacks = callbacks,
    output = output,
    dynamics_only = true,
)

initialize!(model.geometry, model)
set!(
    model;
    orography = wave1_orography,
)

simulation = initialize!(model)

run!(
    simulation;
    period = run_period,
    output = true,
)

sigma_full = model.geometry.σ_levels_full
approx_pressure_hpa = sigma_full .* 1000

println("HS94 NH-winter wave-1 orography run finished")
println("truncation = T", truncation)
println("nlayers = ", nlayers)
println("run_period = ", run_period)
println("time_step_at_T31 = ", time_step_at_T31)
println("actual_time_step_seconds = ", model.time_stepping.Δt_sec)
println("forcing = NHWinterHeldSuarez")
println("gamma_km = ", gamma_km)
println("tropopause_pressure_hpa = ", tropopause_pressure_hpa)
println("vortex_edge_latitude = ", vortex_edge_latitude, "N")
println("vortex_width_latitude = ", vortex_width_latitude, " degrees")
println("top_full_level_sigma = ", sigma_full[1])
println("top_full_level_pressure_hpa_if_ps_1000 = ", approx_pressure_hpa[1])
println("vertical_coordinates = custom nonuniform sigma")
println("orography_amplitude_m = ", orography_amplitude_m)
println("zonal_wavenumber = ", zonal_wavenumber)
println("orography_lat_range = ", lat_south, "N to ", lat_north, "N")
println("orography_min_m = ", minimum(model.orography.orography))
println("orography_max_m = ", maximum(model.orography.orography))
println("output_dir = ", model.output.run_path)
println("output_file = ", joinpath(model.output.run_path, model.output.filename))
println("nans_detected = ", simulation.model.feedback.nans_detected)
