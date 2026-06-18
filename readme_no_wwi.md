# No-WWI Implementation Notes

This file records the no-WWI changes added on 2026-06-17 for the HS94 paper-style diabatic-heating run.

## What Was Added

- `SpeedyWeather_no_wwi/`
  - Local copy of the installed SpeedyWeather package.
  - The global package under `~/.julia/packages/SpeedyWeather/...` was not modified.
  - `Manifest.toml` was copied from the Julia v1.12 default environment so the local copy can be activated as its own project.

- `run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl`
  - Copied from `run_hs94_nh_winter_paper_diabatic_heating.jl`.
  - Uses the local `SpeedyWeather_no_wwi` project:

    ```julia
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "SpeedyWeather_no_wwi"); io = devnull)
    using SpeedyWeather
    ```

- `SpeedyWeather_no_wwi/src/dynamics/no_wwi.jl`
  - New no-WWI implementation file.
  - Defines:

    ```julia
    NoNoWWI
    NoWWICorrection
    no_wwi_correction!
    ```

## SpeedyWeather Package Changes

- `SpeedyWeather_no_wwi/src/SpeedyWeather.jl`
  - Added:

    ```julia
    include("dynamics/no_wwi.jl")
    ```

- `SpeedyWeather_no_wwi/src/models/primitive_dry.jl`
  - Added a new model component:

    ```julia
    @component no_wwi::NW = NoNoWWI()
    ```

  - This makes the default behavior unchanged unless the run script explicitly passes `NoWWICorrection(...)`.

- `SpeedyWeather_no_wwi/src/dynamics/tendencies.jl`
  - Added the no-WWI hook after `bernoulli_potential!` and before tracer advection:

    ```julia
    no_wwi_correction!(vars, model)
    ```

  - The hook runs after the original dynamical tendencies are computed.

## No-WWI Correction

The implemented correction is:

```math
\dot X_{\mathrm{new}}
=
\dot X_{\mathrm{original}}
+
W_k
\left(
\overline{\dot X_{\mathrm{ee}}}
-
\dot X_{\mathrm{ee}}
\right)
```

with:

```math
X \in \{\zeta, D, T\}
```

The code does not modify:

```math
\dot{\ln p_s}
```

The zonal-mean replacement is implemented in spectral space:

- `m = 0` modes are left unchanged.
- Only non-zonal `m > 0` eddy-eddy tendency modes are subtracted.

This is equivalent to subtracting the eddy-eddy tendency and adding back its zonal mean.

## Fixed-Level Mask

The mask is fixed by model layer index, not by instantaneous pressure.

Defaults:

```julia
nwwi_mode = :upper
nwwi_top_k = 14
nwwi_bottom_k = 18
```

For `nwwi_mode=upper`, weights are:

```math
W_k = 1,\quad k \le 14
```

```math
W_k = \frac{18-k}{18-14},\quad 14 < k < 18
```

```math
W_k = 0,\quad k \ge 18
```

So:

```math
W_{14}=1,\quad W_{15}=0.75,\quad W_{16}=0.5,\quad W_{17}=0.25,\quad W_{18}=0
```

For `nwwi_mode=lower`, the same transition is reversed:

```math
W_k = 0,\quad k \le 14
```

```math
W_k = \frac{k-14}{18-14},\quad 14 < k < 18
```

```math
W_k = 1,\quad k \ge 18
```

For `nwwi_mode=all`:

```math
W_k = 1
```

for every vertical level.

## Eddy Fields and Terms

The no-WWI workspace declares eddy fields under:

```julia
vars.dynamics.no_wwi
```

Main eddy fields:

```math
u', v', \zeta', D', T', w', (\ln p_s)'
```

Main eddy-eddy contributions included:

```math
v'\zeta',\quad -u'\zeta'
```

```math
-RT'\nabla(\ln p_s)'
```

```math
\nabla \cdot (u'T', v'T')
```

```math
w'u',\quad w'v',\quad w'T'
```

```math
\frac12(u'^2+v'^2)
```

The wind tendencies are converted back to prognostic tendencies through curl/divergence:

```math
\dot\zeta_{\mathrm{ee}},\quad \dot D_{\mathrm{ee}}
```

## Run Interface

The no-WWI run script accepts:

```text
nwwi=true/false
nwwi_mode=upper/all/lower
nwwi_top_k=<integer>
nwwi_bottom_k=<integer>
```

Defaults:

```text
nwwi=true
nwwi_mode=upper
nwwi_top_k=14
nwwi_bottom_k=18
```

Example:

```bash
/Users/jiachengye/.juliaup/bin/julia run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 10 q0=6 m=1
```

Disable no-WWI:

```bash
/Users/jiachengye/.juliaup/bin/julia run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 10 q0=6 m=1 nwwi=false
```

Turn on no-WWI for all levels:

```bash
/Users/jiachengye/.juliaup/bin/julia run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 10 q0=6 m=1 nwwi_mode=all
```

Turn on no-WWI for lower levels below the k=14 to k=18 transition:

```bash
/Users/jiachengye/.juliaup/bin/julia run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 10 q0=6 m=1 nwwi_mode=lower
```

The output directory label includes either:

```text
noWWI_upper_k14to18buffer
noWWI_all_k14to18buffer
noWWI_lower_k14to18buffer
```

or:

```text
noWWIoff
```

## Backups Created

Backups were created before editing copied files:

```text
SpeedyWeather_no_wwi/src/SpeedyWeather.jl.bak_pre_no_wwi_20260617
SpeedyWeather_no_wwi/src/models/primitive_dry.jl.bak_pre_no_wwi_20260617
SpeedyWeather_no_wwi/src/dynamics/tendencies.jl.bak_pre_no_wwi_20260617
run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl.bak_pre_no_wwi_20260617
SpeedyWeather_no_wwi/src/dynamics/no_wwi.jl.bak_pre_nwwi_mode_20260617
run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl.bak_pre_nwwi_mode_20260617
readme_no_wwi.md.bak_pre_nwwi_mode_20260617
```

## Verification Performed

Package load test:

```bash
/Users/jiachengye/.juliaup/bin/julia -e 'import Pkg; Pkg.activate(joinpath(pwd(), "SpeedyWeather_no_wwi"); io=devnull); using SpeedyWeather; println(pathof(SpeedyWeather)); println(isdefined(SpeedyWeather, :NoWWICorrection))'
```

Result:

```text
/Users/jiachengye/Desktop/UChicago/Blockings/SSW_PPT/HS94/SpeedyWeather_no_wwi/src/SpeedyWeather.jl
true
```

Smoke tests:

```bash
/Users/jiachengye/.juliaup/bin/julia run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 1 1 q0=0
```

Result:

```text
segment_nans_detected = false
```

```bash
/Users/jiachengye/.juliaup/bin/julia run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 1 1 q0=0 nwwi=false
```

Result:

```text
segment_nans_detected = false
```

Smoke-test output directories:

```text
hs94_nh_winter_gamma4_wave1_paper_diabatic_heating_q00Kday_p0800hPa_pt200hPa_lat45N_sigphi10p026761414789405_noWWI_k14to18buffer_no_orography_T31L40_top0p5hPa_dt10min_segmented_1dseg_test_daily/
hs94_nh_winter_gamma4_wave1_paper_diabatic_heating_q00Kday_p0800hPa_pt200hPa_lat45N_sigphi10p026761414789405_noWWIoff_no_orography_T31L40_top0p5hPa_dt10min_segmented_1dseg_test_daily/
```

## Notes and Limitations

- This implementation is intended for the current CPU T31L40 runs.
- GPU execution has not been verified.
- The pressure tendency is deliberately left unchanged in v1.
- The vertical mask is layer-index based; it does not follow instantaneous pressure changes.
- The original control script `run_hs94_nh_winter_paper_diabatic_heating.jl` was not modified.

## Remote Runtime Notes

These notes record the remote-server run commands and the progress/logging
conventions used after the no-WWI implementation was pushed.

On the remote server, run from the repository root:

```bash
cd /nas/jiachengye/git_projects/hs94-speedyweather
```

To fetch the no-WWI implementation:

```bash
git pull origin main
```

For first-time setup or after dependency changes:

```bash
julia --project=SpeedyWeather_no_wwi -e 'import Pkg; Pkg.instantiate()'
```

Create a log directory:

```bash
mkdir -p logs
```

Ten-year no-WWI runs with paper-style diabatic heating turned on
(`q0=6`, `m=1`) and 8 Julia threads:

```bash
nohup julia -t 8 run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 10 365 q0=6 m=1 nwwi=true nwwi_mode=upper > logs/no_wwi_upper_10yr_t8.log 2>&1 &
```

```bash
nohup julia -t 8 run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 10 365 q0=6 m=1 nwwi=true nwwi_mode=all > logs/no_wwi_all_10yr_t8.log 2>&1 &
```

```bash
nohup julia -t 8 run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 10 365 q0=6 m=1 nwwi=true nwwi_mode=lower > logs/no_wwi_lower_10yr_t8.log 2>&1 &
```

Ten-year paper-heating control run with no-WWI turned off:

```bash
nohup julia -t 8 run_hs94_nh_winter_paper_diabatic_heating_no_wwi.jl daily 10 365 q0=6 m=1 nwwi=false > logs/paper_heating_control_10yr_t8.log 2>&1 &
```

`logs/*.log` files are shell-level `nohup` logs. They record script output such
as arguments, output directory, segment start/end messages, and NaN status.

`year_0001/progress.txt` is written by SpeedyWeather inside each output
directory. It is the better file for checking time-step progress, percentage,
estimated remaining time, and basic wind/temperature diagnostics.

For yearly runs, the output directory suffix is:

```text
segmented_yearly_daily
```

not:

```text
segmented_365dseg_daily
```

Check no-WWI progress:

```bash
tail -40 hs94_nh_winter_gamma4_wave1_paper_diabatic_heating_*noWWI_upper*_segmented_yearly_daily/year_0001/progress.txt
```

```bash
tail -40 hs94_nh_winter_gamma4_wave1_paper_diabatic_heating_*noWWI_all*_segmented_yearly_daily/year_0001/progress.txt
```

```bash
tail -40 hs94_nh_winter_gamma4_wave1_paper_diabatic_heating_*noWWI_lower*_segmented_yearly_daily/year_0001/progress.txt
```

Check control-run progress:

```bash
tail -40 hs94_nh_winter_gamma4_wave1_paper_diabatic_heating_*noWWIoff*_segmented_yearly_daily/year_0001/progress.txt
```

Check shell logs:

```bash
tail -40 logs/no_wwi_upper_10yr_t8.log
```

```bash
tail -40 logs/no_wwi_all_10yr_t8.log
```

```bash
tail -40 logs/no_wwi_lower_10yr_t8.log
```

```bash
tail -40 logs/paper_heating_control_10yr_t8.log
```

Check active Julia jobs:

```bash
ps -fu jiachengye | grep '[j]ulia.*run_hs94_nh_winter_paper_diabatic_heating_no_wwi'
```

A normal finished segment should include:

```text
segment_nans_detected = false
```

The current output directory labels include:

```text
noWWI_upper_k14to18buffer
noWWI_all_k14to18buffer
noWWI_lower_k14to18buffer
noWWIoff
```
