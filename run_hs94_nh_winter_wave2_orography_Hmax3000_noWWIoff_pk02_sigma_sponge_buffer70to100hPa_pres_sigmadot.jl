# Segmented NH-winter high-top HS94 experiment with wave-2 Hmax=3000 m
# orography, no-WWI off, PK02 sigma levels, a top sponge, and diagnostic
# surface-pressure / sigma-vertical-velocity output.
#
# This script is based on the latest paper-heating no-WWI script:
#   run_hs94_nh_winter_paper_diabatic_heating_no_wwi_pk02_sigma_sponge_buffer70to100hPa.jl
#
# Main changes relative to that script:
# 1. Uses NHWinterHeldSuarez instead of NHWinterHeldSuarezPaperHeating.
# 2. Removes the diabatic-heating perturbation entirely.
# 3. Uses ManualOrography with nonnegative wave-1 topography:
#
#      h(lambda, phi) = 0.5 Hmax [1 + cos(m lambda)]
#                       sin(pi * (phi - 30N) / 60deg)
#
#    inside 30N <= phi <= 90N, and 0 elsewhere. The default meridional peak
#    is therefore centered at 60N, and the height range is 0 to Hmax.
#
# Features retained from the latest no-WWI / PK02 script:
# 1. Local SpeedyWeather_no_wwi package activation.
# 2. Optional no-WWI tendency correction.
# 3. Default no-WWI buffer near 70-100 hPa on the PK02 grid:
#    k=22 to k=24, corresponding to about 65.6-95.7 hPa if ps = 1000 hPa.
# 4. PK02 nonuniform sigma interfaces:
#    sigma_i = (i / 44)^5 for i = 5:44, plus the top interface sigma = 0.
# 5. PK02-style top sponge applied to layers with sigma_full <= 0.0005.
# 6. T31L40, T31 time step of 10 minutes, yearly segmented restarts.
# 7. NetCDF output variables: u, v, temp, mslp, z, pres, and sigmadot.
#    Vorticity is removed.
#
# Usage:
#   julia run_hs94_nh_winter_wave2_orography_Hmax3000_noWWIoff_pk02_sigma_sponge_buffer70to100hPa_pres_sigmadot.jl daily 50 365
#   julia run_hs94_nh_winter_wave2_orography_Hmax3000_noWWIoff_pk02_sigma_sponge_buffer70to100hPa_pres_sigmadot.jl daily 2 1
#
# Arguments:
#   1. output frequency: daily or 6hourly
#   2. number of yearly segments to run, default 10
#   3. optional segment length in days, default 365. This is mainly for short
#      smoke tests, e.g. "daily 2 1 h=1000" runs two one-day segments.
#   Optional keyword arguments:
#      h=<meters>, default 3000.0. This is the maximum topographic height.
#      m=<integer>, default 2
#      lat_south=<degrees_north>, default 30
#      lat_north=<degrees_north>, default 90
#      gamma=<K/km>, default 4.0
#      nwwi=<true/false>, default false
#      nwwi_mode=<upper/all/lower>, default upper
#      nwwi_top_k=<integer>, default 22
#      nwwi_bottom_k=<integer>, default 24

import Pkg

Pkg.activate(joinpath(@__DIR__, "SpeedyWeather_no_wwi"); io = devnull)

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

Base.@kwdef mutable struct SigmaDotOutput <: SpeedyWeather.AbstractOutputVariable
    name::String = "sigmadot"
    unit::String = "s^-1"
    long_name::String = "sigma-coordinate vertical velocity d sigma/dt"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

SpeedyWeather.path(::SigmaDotOutput, simulation) = simulation.variables.dynamics.w

struct CombinedDrag{D <: Tuple} <: SpeedyWeather.AbstractDrag
    drags::D
end

CombinedDrag(drags...) = CombinedDrag(drags)

function SpeedyWeather.initialize!(drag::CombinedDrag, model::SpeedyWeather.AbstractModel)
    for component in drag.drags
        SpeedyWeather.initialize!(component, model)
    end

    return nothing
end

function SpeedyWeather.drag!(vars, drag::CombinedDrag, lf::Integer, model::SpeedyWeather.AbstractModel)
    for component in drag.drags
        SpeedyWeather.drag!(vars, component, lf, model)
    end

    return nothing
end

Base.@kwdef struct TopSigmaSpongeDrag{NF, VectorType} <: SpeedyWeather.AbstractDrag
    sigma_sp::NF = 0.0005
    kmax::NF = 0.5 / 86400
    drag_coefs::VectorType
end

function TopSigmaSpongeDrag(spectral_grid::SpectralGrid; kwargs...)
    drag_coefs = zeros(spectral_grid.NF, spectral_grid.nlayers)
    return TopSigmaSpongeDrag{spectral_grid.NF, spectral_grid.VectorType}(;
        drag_coefs = drag_coefs,
        kwargs...,
    )
end

function SpeedyWeather.initialize!(drag::TopSigmaSpongeDrag, model::SpeedyWeather.PrimitiveEquation)
    sigma_levels_full = model.geometry.σ_levels_full
    (; sigma_sp, kmax, drag_coefs) = drag

    @. drag_coefs = ifelse(
        sigma_levels_full <= sigma_sp,
        kmax * ((sigma_sp - sigma_levels_full) / sigma_sp)^2,
        zero(kmax),
    )

    return nothing
end

function SpeedyWeather.drag!(vars, drag::TopSigmaSpongeDrag, lf::Integer, model::SpeedyWeather.AbstractModel)
    (; u, v) = vars.grid
    Fu = vars.tendencies.grid.u
    Fv = vars.tendencies.grid.v

    c = vars.prognostic.scale[]
    @. Fu -= c * drag.drag_coefs' .* u
    @. Fv -= c * drag.drag_coefs' .* v

    return nothing
end

function usage_and_exit()
    println("Usage:")
    println("  julia run_hs94_nh_winter_wave1_orography_no_wwi_pk02_sigma_sponge_buffer70to100hPa.jl daily 10")
    println("  julia run_hs94_nh_winter_wave1_orography_no_wwi_pk02_sigma_sponge_buffer70to100hPa.jl 6hourly 10 h=1000")
    println("  julia run_hs94_nh_winter_wave1_orography_no_wwi_pk02_sigma_sponge_buffer70to100hPa.jl daily 2 1 h=1000 gamma=4")
    println("")
    println("Arguments: output_frequency n_segments [segment_days] [h=<meters>] [m=<integer>] [gamma=<K/km>] [nwwi=true]")
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

function parse_keyword_bool(keyword_args, keys, default, name)
    for key in keys
        if haskey(keyword_args, key)
            value = lowercase(strip(keyword_args[key]))

            value in ("true", "t", "yes", "y", "1", "on") && return true
            value in ("false", "f", "no", "n", "0", "off") && return false

            error("Invalid $(name): $(keyword_args[key]). Use true or false.")
        end
    end

    return default
end

function parse_keyword_symbol(keyword_args, keys, default, allowed, name)
    for key in keys
        if haskey(keyword_args, key)
            value = Symbol(lowercase(strip(keyword_args[key])))

            value in allowed && return value

            error("Invalid $(name): $(keyword_args[key]). Use one of $(collect(allowed)).")
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

gamma_km = parse_keyword_float(keyword_args, ("gamma", "gamma_km"), 4.0, "gamma")
tropopause_pressure_hpa = 100.0
vortex_edge_latitude = 50.0
vortex_width_latitude = 10.0

orography_amplitude_m = parse_keyword_float(
    keyword_args,
    ("h", "height", "amp", "amplitude", "orography_amplitude_m"),
    3000.0,
    "orography maximum height",
)
zonal_wavenumber = parse_keyword_int(keyword_args, ("m", "k", "wavenumber", "zonal_wavenumber"), 2, "zonal wavenumber")
lat_south = parse_keyword_float(keyword_args, ("lat_south", "south", "phis"), 30.0, "southern orography latitude")
lat_north = parse_keyword_float(keyword_args, ("lat_north", "north", "phin"), 90.0, "northern orography latitude")

nwwi_enabled = parse_keyword_bool(keyword_args, ("nwwi", "no_wwi"), false, "nwwi")
nwwi_mode = parse_keyword_symbol(
    keyword_args,
    ("nwwi_mode", "no_wwi_mode"),
    :upper,
    Set([:upper, :all, :lower]),
    "no-WWI mode",
)
default_nwwi_top_k = 22
default_nwwi_bottom_k = 24
nwwi_top_k = parse_keyword_int(keyword_args, ("nwwi_top_k", "no_wwi_top_k"), default_nwwi_top_k, "no-WWI top layer")
nwwi_bottom_k = parse_keyword_int(keyword_args, ("nwwi_bottom_k", "no_wwi_bottom_k"), default_nwwi_bottom_k, "no-WWI bottom layer")

orography_amplitude_m >= 0 || error("orography height h must be nonnegative.")
zonal_wavenumber > 0 || error("zonal wavenumber must be positive.")
lat_north > lat_south || error("lat_north must be greater than lat_south.")
lat_south >= -90 || error("lat_south must be at least -90.")
lat_north <= 90 || error("lat_north must not exceed 90.")
nwwi_top_k >= 0 || error("nwwi_top_k must be nonnegative.")
nwwi_bottom_k > nwwi_top_k || error("nwwi_bottom_k must be greater than nwwi_top_k.")
nwwi_bottom_k <= nlayers || error("nwwi_bottom_k must not exceed nlayers=$(nlayers).")

allowed_keywords = Set([
    "h",
    "height",
    "amp",
    "amplitude",
    "orography_amplitude_m",
    "gamma",
    "gamma_km",
    "m",
    "k",
    "wavenumber",
    "zonal_wavenumber",
    "lat_south",
    "south",
    "phis",
    "lat_north",
    "north",
    "phin",
    "nwwi",
    "no_wwi",
    "nwwi_mode",
    "no_wwi_mode",
    "nwwi_top_k",
    "no_wwi_top_k",
    "nwwi_bottom_k",
    "no_wwi_bottom_k",
])
unknown_keywords = setdiff(Set(keys(keyword_args)), allowed_keywords)
isempty(unknown_keywords) || error("Unknown keyword arguments: $(collect(unknown_keywords))")

gamma_label = number_label(gamma_km)
height_label = number_label(orography_amplitude_m)
lat_south_label = number_label(lat_south)
lat_north_label = number_label(lat_north)
nwwi_label = nwwi_enabled ?
    "noWWI_$(String(nwwi_mode))_p70to100hPa_k$(nwwi_top_k)to$(nwwi_bottom_k)buffer" :
    "noWWIoff"

parent_output_dir = "hs94_nh_winter_gamma$(gamma_label)_wave$(zonal_wavenumber)_nonnegative_orography_Hmax$(height_label)m_lat$(lat_south_label)to$(lat_north_label)N_$(nwwi_label)_pk02sigma_sponge0p5hPa_buffer70to100hPa_T31L40_dt10min_segmented_$(segment_label)_$(output_frequency)_pres_sigmadot"

pk02_sigma_n = 44
pk02_sigma_first_i = 5
sponge_sigma_threshold = 0.0005
sponge_kmax_per_day = 0.5

sigma_half = vcat(
    [0.0],
    [(i / pk02_sigma_n)^5 for i in pk02_sigma_first_i:pk02_sigma_n],
)
length(sigma_half) == nlayers + 1 ||
    error("PK02 sigma_half length must be nlayers + 1.")

orography_center_latitude = 0.5 * (lat_south + lat_north)

function wave1_orography(longitude, latitude)
    if lat_south <= latitude <= lat_north
        meridional_envelope = sind(
            180 * (latitude - lat_south) / (lat_north - lat_south),
        )
        return 0.5 *
            orography_amplitude_m *
            (1 + cosd(zonal_wavenumber * longitude)) *
            meridional_envelope
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
    add!(output, SpeedyWeather.SurfacePressureOutput())
    add!(output, GeopotentialHeightOutput())
    add!(output, SigmaDotOutput())
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
    drag = CombinedDrag(
        LinearDrag(spectral_grid),
        TopSigmaSpongeDrag(
            spectral_grid;
            sigma_sp = sponge_sigma_threshold,
            kmax = sponge_kmax_per_day / 86400,
        ),
    )
    orography = ManualOrography(spectral_grid)
    no_wwi = nwwi_enabled ?
        NoWWICorrection(
            spectral_grid;
            active = true,
            mode = nwwi_mode,
            top_k = nwwi_top_k,
            bottom_k = nwwi_bottom_k,
        ) :
        NoNoWWI()
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
            no_wwi = no_wwi,
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
        no_wwi = no_wwi,
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
sponge_active_layers = findall(sigma -> sigma <= sponge_sigma_threshold, geometry.σ_levels_full)
nwwi_top_pressure_hpa_if_ps_1000 = geometry.σ_levels_full[nwwi_top_k] * 1000
nwwi_bottom_pressure_hpa_if_ps_1000 = geometry.σ_levels_full[nwwi_bottom_k] * 1000

println("HS94 segmented NH-winter wave-1 orography run with no-WWI, PK02 sigma, top sponge, and 70-100 hPa no-WWI buffer")
println("parent_output_dir = ", abspath(parent_output_dir))
println("n_segments = ", n_segments)
println("segment_days = ", segment_days)
println("output_frequency = ", output_frequency)
println("output_interval = ", output_interval)
println("restart_policy = segment N starts from year_(N-1)/restart.jld2")
println("forcing = NHWinterHeldSuarez")
println("diabatic_heating = none")
println("gamma_km = ", gamma_km)
println("tropopause_pressure_hpa = ", tropopause_pressure_hpa)
println("vortex_edge_latitude = ", vortex_edge_latitude, "N")
println("vortex_width_latitude = ", vortex_width_latitude, " degrees")
println("orography = ManualOrography")
println("orography_shape = nonnegative wave cosine")
println("orography_max_height_m = ", orography_amplitude_m)
println("zonal_wavenumber = ", zonal_wavenumber)
println("orography_lat_range = ", lat_south, "N to ", lat_north, "N")
println("orography_center_latitude = ", orography_center_latitude, "N")
println("no_wwi_enabled = ", nwwi_enabled)
println("no_wwi_mode = ", nwwi_mode)
println("no_wwi_top_k = ", nwwi_top_k)
println("no_wwi_bottom_k = ", nwwi_bottom_k)
println("no_wwi_top_pressure_hpa_if_ps_1000 = ", nwwi_top_pressure_hpa_if_ps_1000)
println("no_wwi_bottom_pressure_hpa_if_ps_1000 = ", nwwi_bottom_pressure_hpa_if_ps_1000)
println("no_wwi_mask = upper: W=1 above top_k; lower: W=1 below bottom_k; all: W=1 everywhere")
println("vertical_coordinates = PK02 sigma_i=(i/", pk02_sigma_n, ")^5 for i=", pk02_sigma_first_i, ":", pk02_sigma_n, ", plus top sigma=0")
println("top_sponge_sigma_threshold = ", sponge_sigma_threshold)
println("top_sponge_kmax_per_day = ", sponge_kmax_per_day)
println("top_sponge_active_layers = ", sponge_active_layers)

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
    println("segment_orography_min_m = ", minimum(model.orography.orography))
    println("segment_orography_max_m = ", maximum(model.orography.orography))

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

println("HS94 segmented NH-winter wave-1 orography run with no-WWI, PK02 sigma, top sponge, and 70-100 hPa no-WWI buffer finished")
println("truncation = T", truncation)
println("nlayers = ", nlayers)
println("time_step_at_T31 = ", time_step_at_T31)
println("vertical_coordinates = PK02 sigma_i=(i/", pk02_sigma_n, ")^5 for i=", pk02_sigma_first_i, ":", pk02_sigma_n, ", plus top sigma=0")
println("top_full_level_sigma = ", sigma_full[1])
println("top_full_level_pressure_hpa_if_ps_1000 = ", approx_pressure_hpa[1])
println("top_sponge_sigma_threshold = ", sponge_sigma_threshold)
println("top_sponge_kmax_per_day = ", sponge_kmax_per_day)
println("top_sponge_active_layers = ", sponge_active_layers)
println("gamma_km = ", gamma_km)
println("tropopause_pressure_hpa = ", tropopause_pressure_hpa)
println("orography_max_height_m = ", orography_amplitude_m)
println("zonal_wavenumber = ", zonal_wavenumber)
println("orography_lat_range = ", lat_south, "N to ", lat_north, "N")
println("orography_center_latitude = ", orography_center_latitude, "N")
println("no_wwi_enabled = ", nwwi_enabled)
println("no_wwi_mode = ", nwwi_mode)
println("no_wwi_top_k = ", nwwi_top_k)
println("no_wwi_bottom_k = ", nwwi_bottom_k)
println("no_wwi_top_pressure_hpa_if_ps_1000 = ", nwwi_top_pressure_hpa_if_ps_1000)
println("no_wwi_bottom_pressure_hpa_if_ps_1000 = ", nwwi_bottom_pressure_hpa_if_ps_1000)
println("output_variables = u, v, temp, mslp, z, pres, sigmadot")
println("parent_output_dir = ", abspath(parent_output_dir))
println("nans_status = see segment_nans_detected lines above")
