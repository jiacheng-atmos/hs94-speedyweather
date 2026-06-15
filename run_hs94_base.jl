# Baseline Held-Suarez 1994 style dry primitive-equation smoke test.
#
# This script is the minimal reference case for this folder:
# 1. Uses SpeedyWeather PrimitiveDryModel.
# 2. Uses the default HeldSuarez Newtonian temperature relaxation.
# 3. Uses LinearDrag for the lower-boundary Rayleigh drag.
# 4. Uses NoOrography, so the lower boundary is flat.
# 5. Uses T31 horizontal spectral truncation and L8 default equally spaced
#    sigma levels. With ps about 1000 hPa, the top full level is about 62.5 hPa.
# 6. Runs for 20 days and only reports whether NaNs were detected.
#
# This is not an SSW experiment: it has no imposed winter polar vortex, no
# high-top stratospheric grid, no sponge layer, and no stationary-wave
# topographic forcing. It is kept as the clean HS94-style baseline.

using SpeedyWeather

truncation = 31
nlayers = 8
run_period = Day(20)

spectral_grid = SpectralGrid(
    trunc = truncation,
    nlayers = nlayers,
)

forcing = HeldSuarez(spectral_grid)
drag = LinearDrag(spectral_grid)
orography = NoOrography(spectral_grid)

model = PrimitiveDryModel(
    spectral_grid;
    forcing = forcing,
    drag = drag,
    orography = orography,
    dynamics_only = true,
)

simulation = initialize!(model)

run!(
    simulation;
    period = run_period,
)

println("HS94 base run finished")
println("truncation = T", truncation)
println("nlayers = ", nlayers)
println("run_period = ", run_period)
println("nans_detected = ", simulation.model.feedback.nans_detected)
