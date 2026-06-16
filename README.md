# HS94 SpeedyWeather Experiment Notes

This folder contains idealized Held-Suarez / NH-winter SpeedyWeather runs used
for testing stationary-wave modulation of the stratospheric polar vortex.

## Repository State

Use this workflow for remote runs:

```bash
git pull origin main
```

Then run scripts from the repository root:

```bash
cd /nas/jiachengye/git_projects/hs94-speedyweather
```

Do not commit model output or conda environments. Large output folders,
`*.nc`, backup files, and local environments should stay outside Git history.

## Main Scripts

### `run_hs94_base.jl`

Baseline Held-Suarez smoke test.

- Dry `PrimitiveDryModel`
- `HeldSuarez` forcing
- `LinearDrag`
- `NoOrography`
- `T31L8`
- Low-top default sigma grid
- Intended as a basic package/model sanity check, not an SSW experiment

### `run_hs94_wave1_orography.jl`

Default low-top HS94 background plus idealized stationary wave-1 topography.

- `T31L8`
- Default Held-Suarez thermal relaxation
- Wave-1 orography from `30N` to `90N`
- Orography form:

```text
h(lambda, phi) = H cos(lambda) sin(pi * (phi - 30N) / 60deg)
```

- `H = 1000 m`
- Zero outside `30N` to `90N`
- Daily NetCDF output
- This is a topographic forcing test, not the high-top NH-winter SSW setup

### `nh_winter_hs_forcing.jl`

Custom forcing component used by the NH-winter experiments.

- Keeps the Held-Suarez tropospheric equilibrium temperature below the
  tropopause pressure threshold
- Replaces the stratospheric equilibrium temperature with an NH-winter
  polar-vortex profile
- Based on Kushner and Polvani style stratospheric setup, mirrored to NH winter
- Default key settings:

```text
tropopause_pressure_hpa = 100
gamma_km = 4
vortex_edge_latitude = 50N
vortex_width_latitude = 10 degrees
```

`gamma_km = 4` is the colder / stronger vortex case.

### `run_hs94_nh_winter_wave1_orography.jl`

Monolithic high-top NH-winter run.

- `T31L40`
- Custom nonuniform sigma grid
- Top full level `sigma = 0.0005`, about `0.5 hPa` if `ps = 1000 hPa`
- `NHWinterHeldSuarez`
- `gamma_km = 4`
- Wave-1 topography from `30N` to `90N`, `H = 1000 m`
- Time step at T31: `20 minutes`
- 20-year integration
- Daily NetCDF output
- Progress printed every 30 days

This script writes one long `output.nc`. It is useful but less restart-friendly
than the segmented script below.

### `run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl`

Monolithic high-top NH-winter run with configurable output interval and optional
short test length.

Usage:

```bash
julia -t 8 run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl daily
```

```bash
julia -t 8 run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl 6hourly
```

Short test run:

```bash
julia -t 8 run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl daily 20
```

The optional second argument is a number of days, so `3650` means 10 years:

```bash
julia -t 8 run_hs94_nh_winter_wave1_orography_20yr_custom_output.jl 6hourly 3650
```

Output variables are restricted to:

```text
u, v, temp, mslp, z
```

`vor` is removed. `z` is geopotential height in meters, computed as
geopotential divided by gravity.

### `run_hs94_nh_winter_wave1_orography_segmented.jl`

Preferred long-run script.

This script runs the high-top NH-winter experiment in yearly segments. Each
segment writes a separate `output.nc` and `restart.jld2` inside one parent
folder. The segmented long-run script currently uses a 10-minute T31 timestep
for better stability.

Smoke test:

```bash
julia -t 8 run_hs94_nh_winter_wave1_orography_segmented.jl daily 2 1
```

This means:

```text
daily = daily NetCDF output
2     = run 2 segments
1     = each segment is 1 day
```

Formal 10-year 6-hourly run:

```bash
nohup julia -t 8 run_hs94_nh_winter_wave1_orography_segmented.jl 6hourly 10 > hs94_segmented_10y_6hourly.log 2>&1 &
```

Formal 10-year daily run:

```bash
nohup julia -t 8 run_hs94_nh_winter_wave1_orography_segmented.jl daily 10 > hs94_segmented_10y_daily.log 2>&1 &
```

Output directory for 10-year 6-hourly run:

```text
/nas/jiachengye/git_projects/hs94-speedyweather/hs94_nh_winter_gamma4_wave1_orography_T31L40_top0p5hPa_dt10min_segmented_yearly_6hourly
```

Inside that directory:

```text
year_0001/output.nc
year_0001/restart.jld2
year_0002/output.nc
year_0002/restart.jld2
...
year_0010/output.nc
year_0010/restart.jld2
```

The parent folder name does not include `10y` or `20y`, so a 10-year run can
later be extended to 20 years by running:

```bash
julia -t 8 run_hs94_nh_winter_wave1_orography_segmented.jl 6hourly 20
```

Completed years are detected by the presence of `restart.jld2` and skipped.

### `run_hs94_nh_winter_wave1_orography_segmented_gamma4_h_ramp_0to2000m.jl`

Preferred ramp experiment for testing how gradually increasing tropospheric
stationary-wave forcing modulates the NH-winter polar vortex.

This script is copied from the latest
`run_hs94_nh_winter_wave1_orography_segmented_custom_gamma_h.jl`, but fixes the
experiment design instead of taking `gamma` and `h` from keyword arguments.

- `gamma_km = 4`
- `T31L40`
- Top full level near `0.5 hPa`
- Time step at T31: `10 minutes`
- Same nonnegative wave-1 orography shape and location: `30N` to `90N`
- Yearly restart segments
- NetCDF variables: `u, v, temp, mslp, z`

The prescribed orography half-amplitude schedule is:

```text
year_0001: h = 0 m
year_0002: h = 0 m
year_0003: h = 250 m
year_0004: h = 500 m
year_0005: h = 750 m
year_0006: h = 1000 m
year_0007: h = 1250 m
year_0008: h = 1500 m
year_0009: h = 1750 m
year_0010: h = 2000 m
```

Here `h` is the half-amplitude in:

```text
topography(lambda, phi) = h * (1 + cos(lambda)) * envelope(phi)
```

So the actual topographic height range is `0` to `2h`. At year 10, `h = 2000
m`, so the maximum mountain height is `4000 m`. If the run is extended beyond
year 10, `h` stays capped at `2000 m`.

Smoke test:

```bash
julia -t 8 run_hs94_nh_winter_wave1_orography_segmented_gamma4_h_ramp_0to2000m.jl daily 2 1
```

Formal 10-year 6-hourly run:

```bash
nohup julia -t 8 run_hs94_nh_winter_wave1_orography_segmented_gamma4_h_ramp_0to2000m.jl 6hourly 10 > hs94_gamma4_h_ramp_10y_6hourly.log 2>&1 &
```

Output directory for the 10-year 6-hourly run:

```text
/nas/jiachengye/git_projects/hs94-speedyweather/hs94_nh_winter_gamma4_wave1_nonnegative_h_ramp0to2000m_Hmax0to4000m_orography_T31L40_top0p5hPa_dt10min_segmented_yearly_6hourly
```

### `run_hs94_nh_winter_wave_heating_constant_custom_q.jl`

Constant-amplitude idealized tropospheric wave-heating experiment.

This script is for a `Qmax` sweep, not for a slow ramp. Each individual run has
constant wave-heating amplitude throughout the integration. It keeps the same
NH-winter high-top background but removes topography:

- `gamma_km = 4`
- `T31L40`
- Top full level near `0.5 hPa`
- Time step at T31: `10 minutes`
- `NoOrography`
- Constant wave heating/cooling centered at `50N` and `500 hPa`
- Default horizontal width: `sigma_phi = 12 deg`
- Default vertical width: `sigma_ln_p = 0.35`
- NetCDF variables: `u, v, temp, mslp, z`

The imposed heating is:

```text
Q' = Qmax * cos(k * (lambda - lambda0))
          * exp(-0.5 * ((phi - 50N) / 12deg)^2)
          * exp(-0.5 * (log(p / 500hPa) / 0.35)^2)
```

For the default widths, the `1/e` range is about `33N` to `67N` and about
`305 hPa` to `820 hPa`.

Smoke test:

```bash
julia -t 8 run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 1 1 qmax=0.5
```

Suggested 10-year daily sweep:

```bash
nohup julia -t 8 run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=0.5 > hs94_wave_heating_Q0p5_daily_10y.log 2>&1 &
```

```bash
nohup julia -t 8 run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=1 > hs94_wave_heating_Q1_daily_10y.log 2>&1 &
```

```bash
nohup julia -t 8 run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=2 > hs94_wave_heating_Q2_daily_10y.log 2>&1 &
```

```bash
nohup julia -t 8 run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=5 > hs94_wave_heating_Q5_daily_10y.log 2>&1 &
```

```bash
nohup julia -t 8 run_hs94_nh_winter_wave_heating_constant_custom_q.jl daily 10 qmax=10 > hs94_wave_heating_Q10_daily_10y.log 2>&1 &
```

Use `6hourly` instead of `daily` if finer output is needed:

```bash
julia -t 8 run_hs94_nh_winter_wave_heating_constant_custom_q.jl 6hourly 10 qmax=2
```

## Restart Logic

`output.nc` is diagnostic output. It is not a restart file.

SpeedyWeather writes restart state to:

```text
restart.jld2
```

The segmented script uses:

```text
year_0001/restart.jld2 -> initial condition for year_0002
year_0002/restart.jld2 -> initial condition for year_0003
```

If a long run is interrupted, rerun the same command. The script skips completed
segments and resumes from the first year whose `restart.jld2` is missing.

If a year folder exists but has no `restart.jld2`, inspect or move that
incomplete folder before rerunning.

## Checking Remote Progress

Check the process:

```bash
ps -p PID -o pid,etime,stat,%cpu,%mem,cmd
```

For the current 10-year 6-hourly run started in this session, the log path was:

```text
/nas/jiachengye/git_projects/hs94-speedyweather/hs94_segmented_10y_6hourly.log
```

Follow the log:

```bash
tail -f hs94_segmented_10y_6hourly.log
```

Exit `tail -f` with `Ctrl-C`, not `Ctrl-Z`.

If the `nohup` log only shows:

```text
nohup: ignoring input
```

the Julia process may still be running with stdout buffered. Check the per-year
SpeedyWeather progress file instead:

```bash
tail -40 hs94_nh_winter_gamma4_wave1_orography_T31L40_top0p5hPa_dt10min_segmented_yearly_6hourly/year_0001/progress.txt
```

Check completed years:

```bash
find hs94_nh_winter_gamma4_wave1_orography_T31L40_top0p5hPa_dt10min_segmented_yearly_6hourly -name restart.jld2
```

Do not run `ncdump` on an `output.nc` file while Julia is actively writing it.
It can return a temporary HDF error. Wait until that year finishes and
`restart.jld2` appears.

## Output Variables

For the custom and segmented NH-winter scripts, NetCDF output contains:

```text
u(time, layer, lat, lon)
v(time, layer, lat, lon)
temp(time, layer, lat, lon)
z(time, layer, lat, lon)
mslp(time, lat, lon)
```

`z` is geopotential height in meters.

`soil_layer` may still appear as a coordinate/dimension in the NetCDF header
because SpeedyWeather defines it by default. It is not a large data variable.

## Useful NetCDF Checks

After a year has completed:

```bash
ncdump -h path/to/year_0001/output.nc
```

Print a point time series for the 10th vertical layer near `30N, 30E`:

```bash
python -c "import xarray as xr; ds=xr.open_dataset('output.nc'); p=ds[['u','v','temp','z']].sel(lat=30, lon=30, method='nearest').isel(layer=9); m=ds[['mslp']].sel(lat=30, lon=30, method='nearest'); out=xr.merge([p,m]); print(out.to_dataframe())"
```

`isel(layer=9)` means the 10th layer.

## Performance Notes

The remote machine has many CPUs, but this T31L40 case is small. Julia with
`-t 8` and `-t 16` can have similar runtime because the model does not use all
threads efficiently at this low resolution. Use `-t 8` as the practical default
unless benchmarking shows otherwise.

For 6-hourly output, each year writes about four times as many time records as
daily output. Runtime may increase, but file size and I/O are the larger risks.

## Known Local Caveat

At the time this README was created, the local working tree had an uncommitted
change in `run_hs94_base.jl` changing `Day(20)` to `Day(40)`. That change was
not part of the pushed segmented-run updates unless committed separately.
