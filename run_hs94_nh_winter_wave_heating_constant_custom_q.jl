# Segmented NH-winter high-top HS94 experiment with constant tropospheric wave heating.
#
# This script is for idealized constant-amplitude heating/cooling sweep
# experiments, not for ramp-up experiments. It keeps the NH-winter high-top
# Held-Suarez background but removes topography and adds a zonal wave thermal
# forcing centered around 500 hPa and 50N.
#
# Main modifications relative to the original high-top NH-winter orography run:
# 1. Uses NHWinterHeldSuarezWaveHeating instead of NHWinterHeldSuarez.
# 2. Uses NoOrography, so the lower boundary is flat.
# 3. Adds constant-amplitude wave heating:
#
#      Q' = Qmax * cos(k * (lambda - lambda0))
#                * exp(-0.5 * ((phi - 50N) / 12deg)^2)
#                * exp(-0.5 * (log(p / 500hPa) / 0.35)^2)
#
# 4. Qmax is specified from the command line in K/day, e.g. qmax=0.5.
# 5. T31L40 with custom nonuniform sigma coordinates and top full level near
#    0.5 hPa for ps = 1000 hPa.
# 6. T31 time step of 10 minutes.
# 7. NetCDF output variables: u, v, temp, mslp, and z. Vorticity is removed.
#
# Usage:
#   julia run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=0.5
#   julia run_hs94_nh_winter_wave_heating_constant_custom_q.jl 6hourly 10 qmax=1
#   julia run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 2 1 qmax=2
#   julia run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=5 k=2
#
# Arguments:
#   1. output frequency: daily or 6hourly
#   2. number of yearly segments to run, default 10
#   3. optional segment length in days, default 365. This is mainly for short
#      smoke tests, e.g. "daily 2 1 qmax=1" runs two one-day segments.
#   Optional keyword arguments:
#      qmax=<K/day>, default 1.0
#      k=<integer>, default 1
#      phase=<degrees>, default 0
#      lat0=<degrees_north>, default 50
#      sigphi=<degrees>, default 12
#      pc=<hPa>, default 500
#      siglnp=<nondimensional>, default 0.35

using SpeedyWeather

import Dates

include("nh_winter_hs_wave_heating_forcing.jl")

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
    println("  julia run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=0.5")
    println("  julia run_hs94_nh_winter_wave_heating_constant_custom_q.jl 6hourly 10 qmax=1")
    println("  julia run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 2 1 qmax=2")
    println("  julia run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=5 k=2")
    println("")
    println("Arguments: output_frequency n_segments [segment_days] qmax=<K/day> [k=<integer>]")
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

function split_args(args)
    positional_args = String[]
    keyword_args = Dict{String, String}()

    for arg in args
        if occursin("=", arg)
            key_value = split(arg, "=", limit = 2)
            key = lowercase(strip(key_value[1]))
            value = strip(key_value[2])
            keyword_args[key] = value
        else
            push!(positional_args, arg)
        end
    end

    return positional_args, keyword_args
end

function parse_keyword_float(keyword_args, keys, default, name)
    for key in keys
        if haskey(keyword_args, key)
            value = tryparse(Float64, keyword_args[key])

            if isnothing(value)
                error("Invalid $(name): $(keyword_args[key]). Use a number.")
            end

            return value
        end
    end

    return default
end

function parse_keyword_int(keyword_args, keys, default, name)
    for key in keys
        if haskey(keyword_args, key)
            value = tryparse(Int, keyword_args[key])

            if isnothing(value)
                error("Invalid $(name): $(keyword_args[key]). Use an integer.")
            end

            return value
        end
    end

    return default
end

function number_label(value)
    if isinteger(value)
        return string(Int(round(value)))
    end

    return replace(string(value), "." => "p", "-" => "m")
end

function year_folder_name(segment_index)
    return "year_" * lpad(string(segment_index), 4, '0')
end

function restart_file(parent_output_dir, segment_index)
    return joinpath(parent_output_dir, year_folder_name(segment_index), "restart.jld2")
end

truncation = 31
nlayers = 40
default_segments = 10
default_segment_days = 365
time_step_at_T31 = Dates.Minute(10)

positional_args, keyword_args = split_args(ARGS)

output_frequency, output_interval = parse_output_schedule(positional_args)
n_segments = parse_positive_int(positional_args, 2, default_segments, "number of segments")
segment_days = parse_positive_int(positional_args, 3, default_segment_days, "segment length in days")
segment_period = Dates.Day(segment_days)

segment_label = segment_days == 365 ? "yearly" : "$(segment_days)dseg_test"

gamma_km = 4.0
tropopause_pressure_hpa = 100.0
vortex_edge_latitude = 50.0
vortex_width_latitude = 10.0

qmax_kday = parse_keyword_float(keyword_args, ("qmax", "q", "qmax_kday"), 1.0, "Qmax")
zonal_wavenumber = parse_keyword_int(keyword_args, ("k", "wavenumber", "zonal_wavenumber"), 1, "zonal wavenumber")
heating_longitude_phase_deg = parse_keyword_float(keyword_args, ("phase", "lon0", "lambda0"), 0.0, "longitude phase")
heating_latitude_center_deg = parse_keyword_float(keyword_args, ("lat0", "phi0"), 50.0, "latitude center")
heating_latitude_sigma_deg = parse_keyword_float(keyword_args, ("sigphi", "sigma_phi"), 12.0, "latitude sigma")
heating_pressure_center_hpa = parse_keyword_float(keyword_args, ("pc", "p0", "pressure_center_hpa"), 500.0, "pressure center")
heating_sigma_ln_pressure = parse_keyword_float(keyword_args, ("siglnp", "sigma_lnp"), 0.35, "log-pressure sigma")

qmax_kday >= 0 || error("Qmax must be nonnegative.")
zonal_wavenumber > 0 || error("zonal wavenumber must be positive.")
heating_latitude_sigma_deg > 0 || error("latitude sigma must be positive.")
heating_pressure_center_hpa > 0 || error("pressure center must be positive.")
heating_sigma_ln_pressure > 0 || error("log-pressure sigma must be positive.")

allowed_keywords = Set([
    "qmax",
    "q",
    "qmax_kday",
    "k",
    "wavenumber",
    "zonal_wavenumber",
    "phase",
    "lon0",
    "lambda0",
    "lat0",
    "phi0",
    "sigphi",
    "sigma_phi",
    "pc",
    "p0",
    "pressure_center_hpa",
    "siglnp",
    "sigma_lnp",
])
unknown_keywords = setdiff(Set(keys(keyword_args)), allowed_keywords)
isempty(unknown_keywords) || error("Unknown keyword arguments: $(collect(unknown_keywords))")

gamma_label = number_label(gamma_km)
qmax_label = number_label(qmax_kday)
lat_label = number_label(heating_latitude_center_deg)
sigphi_label = number_label(heating_latitude_sigma_deg)
pc_label = number_label(heating_pressure_center_hpa)
siglnp_label = number_label(heating_sigma_ln_pressure)

parent_output_dir = "hs94_nh_winter_gamma$(gamma_label)_wave$(zonal_wavenumber)_constant_heating_Qmax$(qmax_label)Kday_pc$(pc_label)hPa_lat$(lat_label)N_sigphi$(sigphi_label)_siglnp$(siglnp_label)_no_orography_T31L40_top0p5hPa_dt10min_segmented_$(segment_label)_$(output_frequency)"

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
    forcing = NHWinterHeldSuarezWaveHeating(
        spectral_grid;
        gamma_km = gamma_km,
        tropopause_pressure_hpa = tropopause_pressure_hpa,
        vortex_edge_latitude = vortex_edge_latitude,
        vortex_width_latitude = vortex_width_latitude,
        qmax_kday = qmax_kday,
        zonal_wavenumber = zonal_wavenumber,
        heating_longitude_phase_deg = heating_longitude_phase_deg,
        heating_latitude_center_deg = heating_latitude_center_deg,
        heating_latitude_sigma_deg = heating_latitude_sigma_deg,
        heating_pressure_center_hpa = heating_pressure_center_hpa,
        heating_sigma_ln_pressure = heating_sigma_ln_pressure,
    )
    time_stepping = Leapfrog(
        spectral_grid;
        Δt_at_T31 = time_step_at_T31,
    )
    drag = LinearDrag(spectral_grid)
    orography = NoOrography(spectral_grid)
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

println("HS94 segmented NH-winter constant wave-heating run")
println("parent_output_dir = ", abspath(parent_output_dir))
println("n_segments = ", n_segments)
println("segment_days = ", segment_days)
println("output_frequency = ", output_frequency)
println("output_interval = ", output_interval)
println("restart_policy = segment N starts from year_(N-1)/restart.jld2")
println("gamma_km = ", gamma_km)
println("orography = NoOrography")
println("qmax_kday = ", qmax_kday)
println("zonal_wavenumber = ", zonal_wavenumber)
println("heating_longitude_phase_deg = ", heating_longitude_phase_deg)
println("heating_latitude_center_deg = ", heating_latitude_center_deg)
println("heating_latitude_sigma_deg = ", heating_latitude_sigma_deg)
println("heating_pressure_center_hpa = ", heating_pressure_center_hpa)
println("heating_sigma_ln_pressure = ", heating_sigma_ln_pressure)

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
one_over_e_lat_width_deg = sqrt(2) * heating_latitude_sigma_deg
one_over_e_pressure_low_hpa = heating_pressure_center_hpa * exp(-sqrt(2) * heating_sigma_ln_pressure)
one_over_e_pressure_high_hpa = heating_pressure_center_hpa * exp(sqrt(2) * heating_sigma_ln_pressure)

println("HS94 segmented NH-winter constant wave-heating run finished")
println("truncation = T", truncation)
println("nlayers = ", nlayers)
println("time_step_at_T31 = ", time_step_at_T31)
println("top_full_level_sigma = ", sigma_full[1])
println("top_full_level_pressure_hpa_if_ps_1000 = ", approx_pressure_hpa[1])
println("gamma_km = ", gamma_km)
println("tropopause_pressure_hpa = ", tropopause_pressure_hpa)
println("qmax_kday = ", qmax_kday)
println("zonal_wavenumber = ", zonal_wavenumber)
println("heating_center = ", heating_latitude_center_deg, "N, ", heating_pressure_center_hpa, " hPa")
println("heating_1_over_e_lat_range = ", heating_latitude_center_deg - one_over_e_lat_width_deg, "N to ", heating_latitude_center_deg + one_over_e_lat_width_deg, "N")
println("heating_1_over_e_pressure_range_hpa = ", one_over_e_pressure_low_hpa, " to ", one_over_e_pressure_high_hpa)
println("output_variables = u, v, temp, mslp, z")
println("parent_output_dir = ", abspath(parent_output_dir))
println("nans_status = see segment_nans_detected lines above")
