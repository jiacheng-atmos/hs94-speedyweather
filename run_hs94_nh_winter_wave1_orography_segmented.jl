# Segmented NH-winter high-top HS94 experiment with yearly restart files.
#
# Usage:
#   julia run_hs94_nh_winter_wave1_orography_segmented.jl daily 10
#   julia run_hs94_nh_winter_wave1_orography_segmented.jl 6hourly 20
#   julia run_hs94_nh_winter_wave1_orography_segmented.jl daily 2 1
#
# Arguments:
#   1. output frequency: daily or 6hourly
#   2. number of yearly segments to run, default 20
#   3. optional segment length in days, default 365. This is mainly for short
#      smoke tests, e.g. "daily 2 1" runs two one-day restart segments.
#
# The script creates one parent output folder. Inside it, each segment gets one
# subfolder. For normal 365-day segments, the parent folder name does not
# include the requested number of segments, so a 10-year run can later be
# extended to 20 years in the same folder:
#
#   year_0001/output.nc
#   year_0001/restart.jld2
#   year_0002/output.nc
#   year_0002/restart.jld2
#   ...
#
# Segment 1 starts from the default model initial conditions. Segment N>1 starts
# from year_(N-1)/restart.jld2 using SpeedyWeather.StartFromFile. This means an
# interrupted run can resume from the last completed segment. The diagnostic
# output.nc file is not a restart file; restart.jld2 is the restart state.
#
# Physics and numerics match run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl:
# 1. PrimitiveDryModel with NHWinterHeldSuarez forcing.
# 2. gamma_km = 4.0 K/km.
# 3. T31L40 with custom nonuniform sigma coordinates and top full level near
#    0.5 hPa for ps = 1000 hPa.
# 4. Stationary wave-1 orography from 30N to 90N.
# 5. T31 time step of 20 minutes.
# 6. NetCDF output variables: u, v, temp, mslp, and z. Vorticity is removed.

using SpeedyWeather

import Dates

include("nh_winter_hs_forcing.jl")

Base.@kwdef mutable struct PrintSegmentProgress <: SpeedyWeather.AbstractCallback
    interval_days::Float64 = 30.0
    next_day::Float64 = 0.0
    total_days::Float64 = 0.0
    segment_index::Int = 1
    total_segments::Int = 1
end

function SpeedyWeather.initialize!(
    callback::PrintSegmentProgress,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    callback.next_day = callback.interval_days
    callback.total_days = Dates.Second(vars.prognostic.clock.period).value / 86400
    println(
        "segment ",
        callback.segment_index,
        " / ",
        callback.total_segments,
        ": day 0.0 / ",
        callback.total_days,
    )
    return nothing
end

function SpeedyWeather.callback!(
    callback::PrintSegmentProgress,
    vars::SpeedyWeather.Variables,
    model::SpeedyWeather.AbstractModel,
)
    elapsed_days = Dates.Millisecond(vars.prognostic.clock.time - vars.prognostic.clock.start).value / 86400000

    while elapsed_days + 1.0e-9 >= callback.next_day
        println(
            "segment ",
            callback.segment_index,
            " / ",
            callback.total_segments,
            ": day ",
            callback.next_day,
            " / ",
            callback.total_days,
        )
        callback.next_day += callback.interval_days
    end

    return nothing
end

SpeedyWeather.finalize!(::PrintSegmentProgress, args...) = nothing

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

function usage_and_exit()
    println("Usage:")
    println("  julia run_hs94_nh_winter_wave1_orography_segmented.jl daily 10")
    println("  julia run_hs94_nh_winter_wave1_orography_segmented.jl 6hourly 20")
    println("  julia run_hs94_nh_winter_wave1_orography_segmented.jl daily 2 1")
    println("")
    println("Arguments: output_frequency n_segments [segment_days]")
    exit()
end

function parse_output_schedule(args)
    isempty(args) && return "daily", Dates.Day(1)

    option = lowercase(strip(args[1]))

    option in ("-h", "--help", "help") && usage_and_exit()
    option in ("daily", "day", "1day", "1d") && return "daily", Dates.Day(1)
    option in ("6hourly", "6-hourly", "6h", "6hr", "6hrs", "6hour", "6hours") &&
        return "6hourly", Dates.Hour(6)

    error("Unknown output frequency: $(args[1]). Use daily or 6hourly.")
end

function parse_positive_int(args, index, default, name)
    length(args) < index && return default

    value = tryparse(Int, args[index])

    if isnothing(value) || value <= 0
        error("Invalid $(name): $(args[index]). Use a positive integer.")
    end

    return value
end

function year_folder_name(segment_index)
    return "year_" * lpad(string(segment_index), 4, '0')
end

function restart_file(parent_output_dir, segment_index)
    return joinpath(parent_output_dir, year_folder_name(segment_index), "restart.jld2")
end

truncation = 31
nlayers = 40
default_segments = 20
default_segment_days = 365
time_step_at_T31 = Dates.Minute(20)

output_frequency, output_interval = parse_output_schedule(ARGS)
n_segments = parse_positive_int(ARGS, 2, default_segments, "number of segments")
segment_days = parse_positive_int(ARGS, 3, default_segment_days, "segment length in days")
segment_period = Dates.Day(segment_days)

segment_label = segment_days == 365 ? "yearly" : "$(segment_days)dseg_test"
parent_output_dir = "hs94_nh_winter_gamma4_wave1_orography_T31L40_top0p5hPa_dt20min_segmented_$(segment_label)_$(output_frequency)"

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

function wave1_orography(longitude, latitude)
    if lat_south <= latitude <= lat_north
        meridional_envelope = sind(180 * (latitude - lat_south) / (lat_north - lat_south))
        return orography_amplitude_m * cosd(zonal_wavenumber * longitude) * meridional_envelope
    end

    return 0.0
end

function make_output(spectral_grid, segment_index)
    output = NetCDFOutput(
        spectral_grid,
        PrimitiveDry;
        path = parent_output_dir,
        run_prefix = "year",
        run_number = segment_index,
        run_digits = 4,
        interval = output_interval,
    )
    delete!(output, :vor)
    add!(output, GeopotentialHeightOutput())
    return output
end

function make_model(spectral_grid, geometry, segment_index)
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
    callbacks = CallbackDict(
        :progress => PrintSegmentProgress(
            interval_days = 30.0,
            segment_index = segment_index,
            total_segments = n_segments,
        ),
    )
    output = make_output(spectral_grid, segment_index)

    if segment_index == 1
        return PrimitiveDryModel(
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
    end

    initial_conditions = (;
        restart = StartFromFile(
            spectral_grid;
            path = parent_output_dir,
            run_prefix = "year",
            run_number = segment_index - 1,
            run_digits = 4,
        ),
    )

    return PrimitiveDryModel(
        spectral_grid;
        geometry = geometry,
        forcing = forcing,
        drag = drag,
        orography = orography,
        time_stepping = time_stepping,
        callbacks = callbacks,
        output = output,
        initial_conditions = initial_conditions,
        dynamics_only = true,
    )
end

mkpath(parent_output_dir)

spectral_grid = SpectralGrid(
    trunc = truncation,
    nlayers = nlayers,
)
vertical_coordinates = SigmaCoordinates(sigma_half)
geometry = Geometry(
    spectral_grid;
    vertical_coordinates = vertical_coordinates,
)

println("HS94 segmented NH-winter wave-1 orography run")
println("parent_output_dir = ", abspath(parent_output_dir))
println("n_segments = ", n_segments)
println("segment_days = ", segment_days)
println("output_frequency = ", output_frequency)
println("output_interval = ", output_interval)
println("restart_policy = segment N starts from year_(N-1)/restart.jld2")

for segment_index in 1:n_segments
    folder = joinpath(parent_output_dir, year_folder_name(segment_index))
    current_restart = restart_file(parent_output_dir, segment_index)

    if isfile(current_restart)
        println("segment ", segment_index, " already complete; found ", current_restart)
        continue
    end

    if isdir(folder)
        error(
            "Segment folder exists without restart.jld2: $(folder). " *
            "Inspect or move this incomplete folder before rerunning.",
        )
    end

    if segment_index > 1
        previous_restart = restart_file(parent_output_dir, segment_index - 1)
        isfile(previous_restart) ||
            error("Cannot start segment $(segment_index); missing previous restart file: $(previous_restart)")
    end

    println("starting segment ", segment_index, " / ", n_segments)

    model = make_model(spectral_grid, geometry, segment_index)
    initialize!(model.geometry, model)
    set!(
        model;
        orography = wave1_orography,
    )

    simulation = initialize!(model)

    run!(
        simulation;
        period = segment_period,
        output = true,
    )

    println("finished segment ", segment_index, " / ", n_segments)
    println("segment_output_dir = ", model.output.run_path)
    println("segment_output_file = ", joinpath(model.output.run_path, model.output.filename))
    println("segment_restart_file = ", joinpath(model.output.run_path, "restart.jld2"))
    println("segment_nans_detected = ", simulation.model.feedback.nans_detected)
end

sigma_full = geometry.σ_levels_full
approx_pressure_hpa = sigma_full .* 1000

println("HS94 segmented NH-winter wave-1 orography run finished")
println("truncation = T", truncation)
println("nlayers = ", nlayers)
println("time_step_at_T31 = ", time_step_at_T31)
println("top_full_level_sigma = ", sigma_full[1])
println("top_full_level_pressure_hpa_if_ps_1000 = ", approx_pressure_hpa[1])
println("gamma_km = ", gamma_km)
println("tropopause_pressure_hpa = ", tropopause_pressure_hpa)
println("orography_amplitude_m = ", orography_amplitude_m)
println("zonal_wavenumber = ", zonal_wavenumber)
println("orography_lat_range = ", lat_south, "N to ", lat_north, "N")
println("output_variables = u, v, temp, mslp, z")
println("parent_output_dir = ", abspath(parent_output_dir))
println("nans_status = see segment_nans_detected lines above")
