# This script is the 20-year NH-winter high-top HS94 experiment with a
# command-line selectable NetCDF output interval.
#
# Usage:
#   julia run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl daily
#   julia run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl 6hourly
#   julia run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl daily 20
#   julia run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl 6hourly 20
#
# If no argument is given, daily output is used for the full 20 years. The
# optional second argument is a test-run length in days.
#
# Relative to the original local HS94 baseline, this experiment adds:
# 1. A dry primitive-equation SpeedyWeather model with HS-style Newtonian
#    temperature relaxation and linear lower-boundary drag.
# 2. NHWinterHeldSuarez thermal forcing: the troposphere keeps the HS profile,
#    while p < 100 hPa uses a perpetual NH-winter polar-vortex equilibrium
#    temperature based on Kushner and Polvani (2004), mirrored to NH winter.
# 3. gamma_km = 4.0 K/km, the colder / stronger-vortex case in Kushner and
#    Polvani's notation.
# 4. A custom nonuniform L40 SigmaCoordinates grid with top full level
#    sigma = 0.0005, about 0.5 hPa if ps = 1000 hPa.
# 5. Stationary wave-1 orography from 30N to 90N with a sine meridional
#    envelope that is zero at 30N and 90N and peaks near 60N.
# 6. A 20 minute T31 time step, used for stability with the high-top grid.
# 7. A 20 year integration by default, with an optional command-line test-run
#    length in days. The output folder name records whether the NetCDF output
#    is daily or 6hourly and whether the run is a shorter test.
# 8. NetCDF output is limited to u, v, temp, mslp, and geopotential height z.
#    The default vorticity output is removed to reduce file size.
#
# This script does not implement restart/checkpointing. The NetCDF output file
# is diagnostic output, not a complete model state for resuming integration.

using SpeedyWeather

import Dates

include("nh_winter_hs_forcing.jl")

Base.@kwdef mutable struct PrintDayProgress <: SpeedyWeather.AbstractCallback
    interval_days::Float64 = 30.0
    next_day::Float64 = 0.0
    total_days::Float64 = 0.0
end

function SpeedyWeather.initialize!(
    callback::PrintDayProgress,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    callback.next_day = callback.interval_days
    callback.total_days = Dates.Second(vars.prognostic.clock.period).value / 86400
    println("day 0.0 / ", callback.total_days)
    return nothing
end

function SpeedyWeather.callback!(
    callback::PrintDayProgress,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    elapsed_days = Dates.Millisecond(vars.prognostic.clock.time - vars.prognostic.clock.start).value / 86400000

    while elapsed_days + 1.0e-9 >= callback.next_day
        println("day ", callback.next_day, " / ", callback.total_days)
        callback.next_day += callback.interval_days
    end

    return nothing
end

SpeedyWeather.finalize!(::PrintDayProgress, args...) = nothing

Base.@kwdef mutable struct GeopotentialHeightOutput <: SpeedyWeather.AbstractOutputVariable
    name::String = "z"
    unit::String = "m"
    long_name::String = "geopotential height"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

function SpeedyWeather.path(::GeopotentialHeightOutput, simulation)
    geopotential_height = simulation.variables.scratch.grid.a
    geopotential_height .= simulation.variables.grid.geopotential
    geopotential_height ./= simulation.model.planet.gravity
    return geopotential_height
end

function parse_output_schedule(args)
    if isempty(args)
        return "daily", Dates.Day(1)
    end

    option = lowercase(strip(args[1]))

    if option in ("-h", "--help", "help")
        println("Usage:")
        println("  julia run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl daily")
        println("  julia run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl 6hourly")
        println("  julia run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl daily 20")
        println("  julia run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl 6hourly 20")
        println("")
        println("The optional second argument is a test-run length in days.")
        exit()
    elseif option in ("daily", "day", "1day", "1d")
        return "daily", Dates.Day(1)
    elseif option in ("6hourly", "6-hourly", "6h", "6hr", "6hrs", "6hour", "6hours")
        return "6hourly", Dates.Hour(6)
    end

    error("Unknown output frequency: $(args[1]). Use daily or 6hourly.")
end

function parse_run_period(args, default_years)
    if length(args) < 2
        return default_years, Dates.Day(365 * default_years), "20y"
    end

    test_days = tryparse(Int, args[2])

    if isnothing(test_days) || test_days <= 0
        error("Invalid test-run length: $(args[2]). Use a positive integer number of days.")
    end

    return default_years, Dates.Day(test_days), "$(test_days)d_test"
end

truncation = 31
nlayers = 40
run_years = 20
time_step_at_T31 = Dates.Minute(20)
output_frequency, output_interval = parse_output_schedule(ARGS)
run_years, run_period, run_length_label = parse_run_period(ARGS, run_years)

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

output_dir_id = "hs94_nh_winter_gamma4_wave1_orography_T31L40_top0p5hPa_dt20min_$(run_length_label)_$(output_frequency)"

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
callbacks = CallbackDict(:progress => PrintDayProgress(interval_days = 30.0))
output = NetCDFOutput(
    spectral_grid,
    PrimitiveDry;
    id = output_dir_id,
    interval = output_interval,
)
delete!(output, :vor)
add!(output, GeopotentialHeightOutput())

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
println("run_years = ", run_years)
println("run_period = ", run_period)
println("run_length_label = ", run_length_label)
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
println("output_frequency = ", output_frequency)
println("output_interval = ", output_interval)
println("output_variables = ", collect(keys(model.output.variables)))
println("orography_min_m = ", minimum(model.orography.orography))
println("orography_max_m = ", maximum(model.orography.orography))
println("output_dir = ", model.output.run_path)
println("output_file = ", joinpath(model.output.run_path, model.output.filename))
println("nans_detected = ", simulation.model.feedback.nans_detected)
