# Held-Suarez 1994 style dry primitive-equation run with idealized wave-1
# topography.
#
# Relative to run_hs94_base.jl, this script keeps the same low-top HS94-style
# background but adds the following modifications:
# 1. Uses the default HeldSuarez Newtonian temperature relaxation, unchanged.
# 2. Uses LinearDrag for the lower-boundary Rayleigh drag, unchanged.
# 3. Uses T31 horizontal truncation and L8 default equally spaced sigma levels.
#    With ps about 1000 hPa, the top full level is about 62.5 hPa.
# 4. Replaces NoOrography with ManualOrography.
# 5. Adds stationary wave-1 orography from 30N to 90N:
#       h(lambda, phi) = H cos(lambda) sin(pi * (phi - 30N) / 60deg)
#    inside the latitude band, and zero elsewhere.
# 6. Uses H = 1000 m, so the mountain peaks near 60N and is zero at 30N and
#    90N.
# 7. Adds a daily progress callback and daily NetCDF output.
#
# This is a topographic stationary-wave forcing test on the default HS94
# background. It is not yet the NH-winter high-top SSW setup; use
# run_hs94_nh_winter_wave1_orography.jl for that experiment.

using SpeedyWeather

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
nlayers = 8
run_period = Day(20)

orography_amplitude_m = 1000.0
zonal_wavenumber = 1
lat_south = 30.0
lat_north = 90.0
output_dir_id = "hs94_wave1_orography_T31L8_20d"

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

forcing = HeldSuarez(spectral_grid)
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
    forcing = forcing,
    drag = drag,
    orography = orography,
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

println("HS94 wave-1 orography run finished")
println("truncation = T", truncation)
println("nlayers = ", nlayers)
println("run_period = ", run_period)
println("orography_amplitude_m = ", orography_amplitude_m)
println("zonal_wavenumber = ", zonal_wavenumber)
println("orography_lat_range = ", lat_south, "N to ", lat_north, "N")
println("orography_min_m = ", minimum(model.orography.orography))
println("orography_max_m = ", maximum(model.orography.orography))
println("output_dir = ", model.output.run_path)
println("output_file = ", joinpath(model.output.run_path, model.output.filename))
println("nans_detected = ", simulation.model.feedback.nans_detected)
