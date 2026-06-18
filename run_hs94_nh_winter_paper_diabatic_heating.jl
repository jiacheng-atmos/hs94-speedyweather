# Segmented NH-winter high-top HS94 experiment with paper-style tropospheric
# diabatic wave heating.
#
# This script keeps the NH-winter high-top Held-Suarez background, removes
# topography, and adds the Lindgren-style tropospheric diabatic heating
# perturbation.
#
# Main modifications relative to the original high-top NH-winter orography run:
# 1. Uses NHWinterHeldSuarezPaperHeating instead of NHWinterHeldSuarez.
# 2. Uses NoOrography, so the lower boundary is flat.
# 3. Adds the paper-style diabatic heating:
#
#      Q' = q0 * sin(m * lambda)
#               * exp(-0.5 * ((phi - 45N) / sigma_phi)^2)
#               * sin(pi * (-log(p / p0)) / log(pt / p0))
#
#    for pt <= p <= p0, and 0 otherwise.
#
# 4. q0 is specified from the command line in K/day, default q0=6.
# 5. T31L40 with custom nonuniform sigma coordinates and top full level near
#    0.5 hPa for ps = 1000 hPa.
# 6. T31 time step of 10 minutes.
# 7. NetCDF output variables: u, v, temp, mslp, and z. Vorticity is removed.
#
# Usage:
#   julia run_hs94_nh_winter_paper_diabatic_heating.jl daily 10
#   julia run_hs94_nh_winter_paper_diabatic_heating.jl 6hourly 10 m=2
#   julia run_hs94_nh_winter_paper_diabatic_heating.jl daily 2 1 q0=6
#
# Arguments:
#   1. output frequency: daily or 6hourly
#   2. number of yearly segments to run, default 10
#   3. optional segment length in days, default 365. This is mainly for short
#      smoke tests, e.g. "daily 2 1 q0=6" runs two one-day segments.
#   Optional keyword arguments:
#      q0=<K/day>, default 6.0
#      m=<integer>, default 1
#      phase=<degrees>, default 0
#      lat0=<degrees_north>, default 45
#      sigphi=<degrees>, default 0.175 radians = about 10.03 degrees
#      p0=<hPa>, default 800
#      pt=<hPa>, default 200

using SpeedyWeather

import Dates

include("nh_winter_hs_paper_heating_forcing.jl")

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
    println("  julia run_hs94_nh_winter_paper_diabatic_heating.jl daily 10")
    println("  julia run_hs94_nh_winter_paper_diabatic_heating.jl 6hourly 10 m=2")
    println("  julia run_hs94_nh_winter_paper_diabatic_heating.jl daily 2 1 q0=6")
    println("")
    println("Arguments: output_frequency n_segments [segment_days] q0=<K/day> [m=<integer>]")
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

paper_sigma_phi_deg = 0.175 * 360 / (2 * pi)

q0_kday = parse_keyword_float(keyword_args, ("q0", "q"), 6.0, "q0")
zonal_wavenumber = parse_keyword_int(keyword_args, ("m", "k", "wavenumber", "zonal_wavenumber"), 1, "zonal wavenumber")
heating_longitude_phase_deg = parse_keyword_float(keyword_args, ("phase", "lon0", "lambda0"), 0.0, "longitude phase")
heating_latitude_center_deg = parse_keyword_float(keyword_args, ("lat0", "phi0"), 45.0, "latitude center")
heating_latitude_sigma_deg = parse_keyword_float(keyword_args, ("sigphi", "sigma_phi"), paper_sigma_phi_deg, "latitude sigma")
heating_pressure_bottom_hpa = parse_keyword_float(keyword_args, ("p0", "pbottom", "pressure_bottom_hpa"), 800.0, "bottom pressure")
heating_pressure_top_hpa = parse_keyword_float(keyword_args, ("pt", "ptop", "pressure_top_hpa"), 200.0, "top pressure")

q0_kday >= 0 || error("q0 must be nonnegative.")
zonal_wavenumber > 0 || error("zonal wavenumber must be positive.")
heating_latitude_sigma_deg > 0 || error("latitude sigma must be positive.")
heating_pressure_bottom_hpa > 0 || error("bottom pressure must be positive.")
heating_pressure_top_hpa > 0 || error("top pressure must be positive.")
heating_pressure_top_hpa < heating_pressure_bottom_hpa ||
    error("top pressure pt must be less than bottom pressure p0.")

allowed_keywords = Set([
    "q",
    "q0",
    "m",
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
    "p0",
    "pbottom",
    "pressure_bottom_hpa",
    "pt",
    "ptop",
    "pressure_top_hpa",
])
unknown_keywords = setdiff(Set(keys(keyword_args)), allowed_keywords)
isempty(unknown_keywords) || error("Unknown keyword arguments: $(collect(unknown_keywords))")

gamma_label = number_label(gamma_km)
q0_label = number_label(q0_kday)
lat_label = number_label(heating_latitude_center_deg)
sigphi_label = number_label(heating_latitude_sigma_deg)
p0_label = number_label(heating_pressure_bottom_hpa)
pt_label = number_label(heating_pressure_top_hpa)

parent_output_dir = "hs94_nh_winter_gamma$(gamma_label)_wave$(zonal_wavenumber)_paper_diabatic_heating_q0$(q0_label)Kday_p0$(p0_label)hPa_pt$(pt_label)hPa_lat$(lat_label)N_sigphi$(sigphi_label)_no_orography_T31L40_top0p5hPa_dt10min_segmented_$(segment_label)_$(output_frequency)"

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
    forcing = NHWinterHeldSuarezPaperHeating(
        spectral_grid;
        gamma_km = gamma_km,
        tropopause_pressure_hpa = tropopause_pressure_hpa,
        vortex_edge_latitude = vortex_edge_latitude,
        vortex_width_latitude = vortex_width_latitude,
        q0_kday = q0_kday,
        zonal_wavenumber = zonal_wavenumber,
        heating_longitude_phase_deg = heating_longitude_phase_deg,
        heating_latitude_center_deg = heating_latitude_center_deg,
        heating_latitude_sigma_deg = heating_latitude_sigma_deg,
        heating_pressure_bottom_hpa = heating_pressure_bottom_hpa,
        heating_pressure_top_hpa = heating_pressure_top_hpa,
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

println("HS94 segmented NH-winter paper diabatic-heating run")
println("parent_output_dir = ", abspath(parent_output_dir))
println("n_segments = ", n_segments)
println("segment_days = ", segment_days)
println("output_frequency = ", output_frequency)
println("output_interval = ", output_interval)
println("restart_policy = segment N starts from year_(N-1)/restart.jld2")
println("gamma_km = ", gamma_km)
println("orography = NoOrography")
println("q0_kday = ", q0_kday)
println("zonal_wavenumber = ", zonal_wavenumber)
println("heating_longitude_phase_deg = ", heating_longitude_phase_deg)
println("heating_latitude_center_deg = ", heating_latitude_center_deg)
println("heating_latitude_sigma_deg = ", heating_latitude_sigma_deg)
println("heating_pressure_bottom_hpa = ", heating_pressure_bottom_hpa)
println("heating_pressure_top_hpa = ", heating_pressure_top_hpa)

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
heating_pressure_peak_magnitude_hpa = sqrt(
    heating_pressure_bottom_hpa * heating_pressure_top_hpa,
)

println("HS94 segmented NH-winter paper diabatic-heating run finished")
println("truncation = T", truncation)
println("nlayers = ", nlayers)
println("time_step_at_T31 = ", time_step_at_T31)
println("top_full_level_sigma = ", sigma_full[1])
println("top_full_level_pressure_hpa_if_ps_1000 = ", approx_pressure_hpa[1])
println("gamma_km = ", gamma_km)
println("tropopause_pressure_hpa = ", tropopause_pressure_hpa)
println("q0_kday = ", q0_kday)
println("zonal_wavenumber = ", zonal_wavenumber)
println("heating_latitude_center = ", heating_latitude_center_deg, "N")
println("heating_1_over_e_lat_range = ", heating_latitude_center_deg - one_over_e_lat_width_deg, "N to ", heating_latitude_center_deg + one_over_e_lat_width_deg, "N")
println("heating_pressure_range_hpa = ", heating_pressure_top_hpa, " to ", heating_pressure_bottom_hpa)
println("heating_pressure_peak_magnitude_hpa = ", heating_pressure_peak_magnitude_hpa)
println("output_variables = u, v, temp, mslp, z")
println("parent_output_dir = ", abspath(parent_output_dir))
println("nans_status = see segment_nans_detected lines above")
