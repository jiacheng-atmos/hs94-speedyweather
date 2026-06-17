abstract type AbstractNoWWI <: AbstractModelComponent end

export NoNoWWI, NoWWICorrection

struct NoNoWWI <: AbstractNoWWI end

Base.@kwdef struct NoWWICorrection <: AbstractNoWWI
    active::Bool = true
    mode::Symbol = :upper
    top_k::Int = 14
    bottom_k::Int = 18
end

NoWWICorrection(::SpectralGrid; kwargs...) = NoWWICorrection(; kwargs...)

Adapt.@adapt_structure NoWWICorrection

variables(::NoNoWWI) = ()

function variables(::NoWWICorrection)
    ns = :no_wwi
    return (
        DynamicsVariable(:u_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:v_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:vorticity_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:divergence_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:temperature_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:w_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:pressure_eddy, Grid2D(), namespace = ns),
        DynamicsVariable(:dpres_dx_eddy, Grid2D(), namespace = ns),
        DynamicsVariable(:dpres_dy_eddy, Grid2D(), namespace = ns),
        DynamicsVariable(:pressure_eddy_spec, Spectral2D(), namespace = ns),
        DynamicsVariable(:dpres_dx_eddy_spec, Spectral2D(), namespace = ns),
        DynamicsVariable(:dpres_dy_eddy_spec, Spectral2D(), namespace = ns),
        DynamicsVariable(:pres_flux_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:div_sum_above_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:pres_flux_sum_above_eddy, Grid3D(), namespace = ns),
        DynamicsVariable(:u_tend_ee, Grid3D(), namespace = ns),
        DynamicsVariable(:v_tend_ee, Grid3D(), namespace = ns),
        DynamicsVariable(:temperature_tend_ee_grid, Grid3D(), namespace = ns),
        DynamicsVariable(:bernoulli_ee_grid, Grid3D(), namespace = ns),
        DynamicsVariable(:vorticity_tend_ee, Spectral3D(), namespace = ns),
        DynamicsVariable(:divergence_tend_ee, Spectral3D(), namespace = ns),
        DynamicsVariable(:temperature_tend_ee, Spectral3D(), namespace = ns),
        DynamicsVariable(:bernoulli_ee, Spectral3D(), namespace = ns),
    )
end

initialize!(::NoNoWWI, ::AbstractModel) = nothing

function initialize!(scheme::NoWWICorrection, model::PrimitiveEquation)
    nlayers = model.geometry.nlayers
    scheme.mode in (:upper, :all, :lower) ||
        error("NoWWICorrection mode must be :upper, :all, or :lower.")
    scheme.top_k >= 0 || error("NoWWICorrection top_k must be nonnegative.")
    scheme.bottom_k > scheme.top_k || error("NoWWICorrection bottom_k must be greater than top_k.")
    scheme.bottom_k <= nlayers || error("NoWWICorrection bottom_k exceeds nlayers=$(nlayers).")
    return nothing
end

no_wwi_correction!(vars::Variables, model::AbstractModel) = nothing
no_wwi_correction!(vars::Variables, model::PrimitiveDry) =
    no_wwi_correction!(vars, model.no_wwi, model)

no_wwi_correction!(vars::Variables, ::NoNoWWI, model::PrimitiveDry) = nothing

function no_wwi_correction!(
        vars::Variables,
        scheme::NoWWICorrection,
        model::PrimitiveDry,
    )
    scheme.active || return nothing

    compute_no_wwi_eddy_fields!(vars, model)
    compute_no_wwi_pressure_gradient!(vars, model)
    compute_no_wwi_pressure_flux_sums!(vars, model)
    compute_no_wwi_vordiv_tendencies!(vars, model)
    compute_no_wwi_temperature_tendency!(vars, model)

    apply_no_wwi_spectral_correction!(
        vars.tendencies.vorticity,
        vars.dynamics.no_wwi.vorticity_tend_ee,
        scheme,
    )
    apply_no_wwi_spectral_correction!(
        vars.tendencies.divergence,
        vars.dynamics.no_wwi.divergence_tend_ee,
        scheme,
    )
    apply_no_wwi_spectral_correction!(
        vars.tendencies.temperature,
        vars.dynamics.no_wwi.temperature_tend_ee,
        scheme,
    )

    return nothing
end

@inline function no_wwi_weight(scheme::NoWWICorrection, k::Integer)
    scheme.mode === :all && return one(Float64)

    if scheme.mode === :upper
        k <= scheme.top_k && return one(Float64)
        k >= scheme.bottom_k && return zero(Float64)
        return (scheme.bottom_k - k) / (scheme.bottom_k - scheme.top_k)
    end

    k <= scheme.top_k && return zero(Float64)
    k >= scheme.bottom_k && return one(Float64)
    return (k - scheme.top_k) / (scheme.bottom_k - scheme.top_k)
end

function apply_no_wwi_spectral_correction!(
        tendency::LowerTriangularArray,
        eddy_eddy_tendency::LowerTriangularArray,
        scheme::NoWWICorrection,
    )
    nlayers = size(tendency, 2)
    for (m1, lms) in enumerate(eachorder(tendency))
        m1 == 1 && continue
        for lm in lms
            for k in 1:nlayers
                weight = no_wwi_weight(scheme, k)
                weight == 0 && continue
                tendency[lm, k] -= weight * eddy_eddy_tendency[lm, k]
            end
        end
    end
    return nothing
end

function compute_no_wwi_eddy_fields!(vars::Variables, model::PrimitiveDry)
    workspace = vars.dynamics.no_wwi
    zonal_eddy!(workspace.u_eddy, vars.grid.u)
    zonal_eddy!(workspace.v_eddy, vars.grid.v)
    zonal_eddy!(workspace.vorticity_eddy, vars.grid.vorticity)
    zonal_eddy!(workspace.divergence_eddy, vars.grid.divergence)
    zonal_eddy!(workspace.temperature_eddy, vars.grid.temperature)
    zonal_eddy!(workspace.w_eddy, vars.dynamics.w)
    zonal_eddy!(workspace.pressure_eddy, vars.grid.pressure)
    return nothing
end

function zonal_eddy!(eddy::AbstractField2D, field::AbstractField2D)
    rings = eachring(field.grid)
    @inbounds for ring in rings
        mean_value = zero(eltype(field))
        for ij in ring
            mean_value += field[ij]
        end
        mean_value /= length(ring)
        for ij in ring
            eddy[ij] = field[ij] - mean_value
        end
    end
    return nothing
end

function zonal_eddy!(eddy::AbstractField3D, field::AbstractField3D)
    rings = eachring(field.grid)
    nlayers = size(field, 2)
    @inbounds for k in 1:nlayers
        for ring in rings
            mean_value = zero(eltype(field))
            for ij in ring
                mean_value += field[ij, k]
            end
            mean_value /= length(ring)
            for ij in ring
                eddy[ij, k] = field[ij, k] - mean_value
            end
        end
    end
    return nothing
end

function compute_no_wwi_pressure_gradient!(vars::Variables, model::PrimitiveDry)
    workspace = vars.dynamics.no_wwi
    S = model.spectral_transform
    scratch_memory = vars.scratch.transform_memory

    pressure_eddy_spec = workspace.pressure_eddy_spec
    dpres_dx_eddy_spec = workspace.dpres_dx_eddy_spec
    dpres_dy_eddy_spec = workspace.dpres_dy_eddy_spec

    transform!(pressure_eddy_spec, workspace.pressure_eddy, scratch_memory, S)
    ∇!(dpres_dx_eddy_spec, dpres_dy_eddy_spec, pressure_eddy_spec, S)
    transform!(workspace.dpres_dx_eddy, dpres_dx_eddy_spec, scratch_memory, S, unscale_coslat = true)
    transform!(workspace.dpres_dy_eddy, dpres_dy_eddy_spec, scratch_memory, S, unscale_coslat = true)

    return nothing
end

function compute_no_wwi_pressure_flux_sums!(vars::Variables, model::PrimitiveDry)
    workspace = vars.dynamics.no_wwi
    (; σ_levels_thick, nlayers) = model.geometry
    (; u_eddy, v_eddy, divergence_eddy, dpres_dx_eddy, dpres_dy_eddy) = workspace
    (; pres_flux_eddy, div_sum_above_eddy, pres_flux_sum_above_eddy) = workspace

    @. pres_flux_eddy = u_eddy * dpres_dx_eddy + v_eddy * dpres_dy_eddy

    @inbounds for ij in eachgridpoint(pres_flux_eddy)
        div_sum = zero(eltype(div_sum_above_eddy))
        pres_flux_sum = zero(eltype(pres_flux_sum_above_eddy))
        for k in 1:nlayers
            div_sum_above_eddy[ij, k] = div_sum
            pres_flux_sum_above_eddy[ij, k] = pres_flux_sum

            Δσk = σ_levels_thick[k]
            div_sum += divergence_eddy[ij, k] * Δσk
            pres_flux_sum += pres_flux_eddy[ij, k] * Δσk
        end
    end

    return nothing
end

function compute_no_wwi_vordiv_tendencies!(vars::Variables, model::PrimitiveDry)
    workspace = vars.dynamics.no_wwi
    (; u_tend_ee, v_tend_ee) = workspace
    fill!(u_tend_ee, 0)
    fill!(v_tend_ee, 0)

    _vertical_advection!(
        u_tend_ee,
        workspace.w_eddy,
        workspace.u_eddy,
        model.geometry.σ_levels_thick,
        model.vertical_advection,
    )
    _vertical_advection!(
        v_tend_ee,
        workspace.w_eddy,
        workspace.v_eddy,
        model.geometry.σ_levels_thick,
        model.vertical_advection,
    )

    add_no_wwi_momentum_terms!(workspace, model)

    S = model.spectral_transform
    scratch_memory = vars.scratch.transform_memory
    u_tend = vars.scratch.a
    v_tend = vars.scratch.b

    transform!(u_tend, u_tend_ee, scratch_memory, S)
    transform!(v_tend, v_tend_ee, scratch_memory, S)

    curl!(workspace.vorticity_tend_ee, u_tend, v_tend, S, add = false)
    divergence!(workspace.divergence_tend_ee, u_tend, v_tend, S, add = false)
    add_no_wwi_bernoulli_tendency!(vars, model)

    return nothing
end

function add_no_wwi_momentum_terms!(workspace, model::PrimitiveDry)
    (; R_dry) = model.atmosphere
    (; coslat⁻¹) = model.geometry
    (; whichring) = workspace.u_eddy.grid

    @inbounds for k in axes(workspace.u_eddy, 2)
        for ij in eachgridpoint(workspace.u_eddy)
            j = whichring[ij]
            temp_pressure_term = R_dry * workspace.temperature_eddy[ij, k]
            workspace.u_tend_ee[ij, k] = (
                workspace.u_tend_ee[ij, k] +
                workspace.v_eddy[ij, k] * workspace.vorticity_eddy[ij, k] -
                temp_pressure_term * workspace.dpres_dx_eddy[ij]
            ) * coslat⁻¹[j]
            workspace.v_tend_ee[ij, k] = (
                workspace.v_tend_ee[ij, k] -
                workspace.u_eddy[ij, k] * workspace.vorticity_eddy[ij, k] -
                temp_pressure_term * workspace.dpres_dy_eddy[ij]
            ) * coslat⁻¹[j]
        end
    end

    return nothing
end

function add_no_wwi_bernoulli_tendency!(vars::Variables, model::PrimitiveDry)
    workspace = vars.dynamics.no_wwi
    half = convert(eltype(workspace.bernoulli_ee_grid), 0.5)
    @. workspace.bernoulli_ee_grid = half * (workspace.u_eddy^2 + workspace.v_eddy^2)

    S = model.spectral_transform
    scratch_memory = vars.scratch.transform_memory
    transform!(workspace.bernoulli_ee, workspace.bernoulli_ee_grid, scratch_memory, S)
    ∇²!(workspace.divergence_tend_ee, workspace.bernoulli_ee, S, add = true, flipsign = true)

    return nothing
end

function compute_no_wwi_temperature_tendency!(vars::Variables, model::PrimitiveDry)
    workspace = vars.dynamics.no_wwi
    temp_tend_grid = workspace.temperature_tend_ee_grid
    fill!(temp_tend_grid, 0)

    _vertical_advection!(
        temp_tend_grid,
        workspace.w_eddy,
        workspace.temperature_eddy,
        model.geometry.σ_levels_thick,
        model.vertical_advection,
    )

    add_no_wwi_temperature_grid_terms!(workspace, model)

    S = model.spectral_transform
    scratch_memory = vars.scratch.transform_memory
    transform!(workspace.temperature_tend_ee, temp_tend_grid, scratch_memory, S)
    add_no_wwi_temperature_flux_divergence!(vars, model)

    return nothing
end

function add_no_wwi_temperature_grid_terms!(workspace, model::PrimitiveDry)
    (; κ) = model.atmosphere
    A = model.adiabatic_conversion.σ_lnp_A
    B = model.adiabatic_conversion.σ_lnp_B

    @inbounds for k in axes(workspace.temperature_eddy, 2)
        Ak = A[k]
        Bk = B[k]
        for ij in eachgridpoint(workspace.temperature_eddy)
            dlnpdt_eddy =
                Ak * (
                    workspace.div_sum_above_eddy[ij, k] +
                    workspace.pres_flux_sum_above_eddy[ij, k]
                ) +
                Bk * (
                    workspace.divergence_eddy[ij, k] +
                    workspace.pres_flux_eddy[ij, k]
                ) +
                workspace.pres_flux_eddy[ij, k]

            workspace.temperature_tend_ee_grid[ij, k] +=
                workspace.temperature_eddy[ij, k] * workspace.divergence_eddy[ij, k] +
                κ * workspace.temperature_eddy[ij, k] * dlnpdt_eddy
        end
    end

    return nothing
end

function add_no_wwi_temperature_flux_divergence!(vars::Variables, model::PrimitiveDry)
    workspace = vars.dynamics.no_wwi
    (; coslat⁻¹) = model.geometry
    (; whichring) = workspace.temperature_eddy.grid

    @inbounds for k in axes(workspace.temperature_eddy, 2)
        for ij in eachgridpoint(workspace.temperature_eddy)
            j = whichring[ij]
            temp_scaled = workspace.temperature_eddy[ij, k] * coslat⁻¹[j]
            workspace.u_tend_ee[ij, k] = workspace.u_eddy[ij, k] * temp_scaled
            workspace.v_tend_ee[ij, k] = workspace.v_eddy[ij, k] * temp_scaled
        end
    end

    S = model.spectral_transform
    scratch_memory = vars.scratch.transform_memory
    u_flux = vars.scratch.a
    v_flux = vars.scratch.b

    transform!(u_flux, workspace.u_tend_ee, scratch_memory, S)
    transform!(v_flux, workspace.v_tend_ee, scratch_memory, S)

    divergence!(workspace.temperature_tend_ee, u_flux, v_flux, S, add = true, flipsign = true)

    return nothing
end
