# CAP Calibration Guide

This guide merges the previous workflow note and data-integration report for the CAP in vitro calibration example.

It documents the current simplified model, the calibration inputs, the make-based workflow, and the current optimization strategy.

## Scope

- **Cell line**: `EGI1`
- **CAP duration**: `5 min`
- **Observed trajectory**: `0 h`, `24 h`, `48 h`, `72 h`
- **Base simulation input**: `input_inVitro.csv`
- **Calibration template**: `input_inVitro_calibration.csv`
- **ABM initial population in calibration template**: `30000`
- **Optimizers**: `differential_evolution`, `dual_annealing`, `nelder_mead`

## Biological interpretation of the current simplified model

The current CAP scenario represents a 3D tumor spheroid exposed to CAP-derived extracellular species:

- `H2O2` as a primary CAP-associated species
- `NO2_` as a secondary CAP-associated species
- `O2` as a baseline metabolic field for growth and survival rules

The current simplified in vitro model intentionally removes the more complicated extracellular dependency logic:

- no migration
- no chemotaxis
- no extracellular `H2O2` / `NO2_` secretion terms for cancer cells
- no RONS-controlled apoptosis, growth, or division dependencies
- CAP species act mainly through intracellular uptake, oxidative stress, and DNA damage

This makes the dominant mechanism:

- extracellular CAP exposure
- intracellular uptake of `H2O2` and `NO2_`
- intracellular oxidative stress accumulation
- DNA damage accumulation and delayed death / growth suppression

## What `input_inVitro.csv` controls

### Global timing and domain

- `number_of_time_steps` and `time_step` set the total simulation duration
- `min_boundary` and `max_boundary` define the diffusion domain
- `diffusion_grid/spatial_resolution` sets the extracellular grid resolution

### Extracellular fields and CAP schedule

- `diffusion_grid/biochemicals` and `substances` declare active fields
- `O2/diffusion_coefficient`, `H2O2/diffusion_coefficient`, and `NO2_/diffusion_coefficient` control extracellular spread
- `H2O2/dissipation_coefficient` and `NO2_/dissipation_coefficient` control extracellular decay
- `CAP/start_step`, `CAP/duration_steps`, `CAP/start_time_h`, and `CAP/duration_h` define the treatment window
- `CAP/H2O2/concentration` and `CAP/NO2_/concentration` define the normalized CAP composition
- `CAP/application_mode` controls how CAP is applied during the active exposure window
- `CAP/post_treatment_mode` controls the CAP-species boundary behavior after the exposure ends

The RONS diffusion coefficients are treated as effective numerical diffusion parameters, not direct aqueous molecular diffusion coefficients.

The current effective values are:

- `H2O2/diffusion_coefficient = 900`
- `NO2_/diffusion_coefficient = 900`

Physical aqueous small-molecule diffusivities can be on the order of `10^-9 m^2/s`, or about `3.6e6 µm^2/h`. Those physical values are intentionally not used in the current explicit diffusion setup because they would require much smaller stable timesteps or a different transport solver.

The CAP transport layer should remain modular:

- current default: explicit BioDynaMo diffusion with effective numerical diffusion
- future option: CAP-specific substepping or an implicit transport solve for physically scaled diffusion

The optional implicit CAP transport path is selected with:

- `CAP/transport_solver = implicit`
- `CAP/implicit_transport/iterations = 50`
- `CAP/implicit_transport/tolerance = 1e-8`
- `CAP/implicit_transport/relaxation = 1.0`

When `CAP/transport_solver = implicit`, BioDynaMo explicit diffusion and dissipation are disabled for `H2O2` and `NO2_`, and the CAP-specific implicit transport step uses the configured species diffusion and dissipation coefficients instead.

For the current example, the CAP concentrations are normalized fractions and should satisfy:

- `0 <= CAP/H2O2/concentration <= 1`
- `0 <= CAP/NO2_/concentration <= 1`
- `CAP/H2O2/concentration + CAP/NO2_/concentration <= 1`

### Cancer-cell behavior block

The cancer-cell block now primarily controls:

- baseline apoptosis, growth, and division
- oxygen-limited baseline behavior
- O2 consumption from the extracellular grid
- intracellular CAP uptake and damage dynamics

The current in-vitro baseline geometry and proliferation setup is:

- `cancer_cell/diameter/min = 10.0`
- `cancer_cell/diameter/max = 16.0`
- `cancer_cell/can_divide/diameter_cutoff = 16.0`
- `cancer_cell/can_grow/probability = 1.0`
- `cancer_cell/can_grow/diameter_rate = 0.20`
- `cancer_cell/can_divide/probability = 0.05`
- `cancer_cell/can_divide/time_window = 150`

This is important because the division cutoff must stay reachable within the configured diameter range.

Baseline proliferation was adjusted using the 15 s CAP condition as a minimum-dose reference, not as an untreated biological control.

## Intracellular CAP mechanism used in the model

At each step for each viable cancer cell:

- local extracellular `H2O2`, `NO2_`, and `O2` are sampled
- `H2O2` and `NO2_` are taken up according to the configured uptake terms
- intracellular oxidative stress is updated
- intracellular DNA damage is updated
- damage can block growth or division and can trigger apoptosis stochastically

The intracellular update is effectively driven by:

$$
\frac{d\,ROS}{dt} = \alpha_{H2O2} \cdot uptake_{H2O2} + \alpha_{NO2} \cdot uptake_{NO2} - k_{scavenge} \cdot antioxidant\_capacity \cdot ROS
$$

and

$$
\frac{d\,Damage}{dt} = k_{induction} \cdot ROS - k_{repair} \cdot Damage
$$

In the ABM calibration workflow, the first-pass fitted parameter set focuses on CAP dose amplitude, intracellular uptake, antioxidant buffering, ROS/RNS weighting, and damage accumulation/repair.

## Which parameters matter most in the current calibration strategy

### Core fitted parameters

These are the parameters currently optimized by the ABM calibration make target:

- `CAP/H2O2/concentration`
- `CAP/NO2_/concentration`
- `cancer_cell/intracellular/uptake/H2O2`
- `cancer_cell/intracellular/uptake/NO2_`
- `cancer_cell/intracellular/antioxidant/k_scavenge`
- `cancer_cell/intracellular/antioxidant/capacity`
- `cancer_cell/intracellular/alpha/H2O2`
- `cancer_cell/intracellular/alpha/NO2_`
- `cancer_cell/intracellular/damage/k_induction`
- `cancer_cell/intracellular/damage/k_repair`
- `cancer_cell/intracellular/damage/threshold`
- `cancer_cell/intracellular/damage/probability`
- `cancer_cell/can_divide/CAP_sensitivity` (β_CAP)

These were chosen because CAP/PAM studies commonly report dose-dependent decreases in cancer-cell survival with ROS/RNS-mediated damage, cell-cycle arrest, and apoptosis.

The CAP-induced proliferation modulation uses a bounded function:
```
p_div(t_CAP) = p_0 × [1 + β_CAP × t_CAP / (K_CAP + t_CAP)]
```
where β_CAP (`cancer_cell/can_divide/CAP_sensitivity`) is the strength of CAP-induced division increase and K_CAP (`cancer_cell/can_divide/CAP_saturation_time`, fixed at 60 s initially) is the saturation time. This captures the experimental observation that longer CAP treatments can increase division probability while preventing overfitting with per-duration parameters.

### Fixed during calibration

All other fields in `input_inVitro_calibration.csv` remain fixed during an optimizer-driven calibration run, including:

- the CAP treatment duration and timing fields generated by `prepare_calibration_input.py`

The fitted CAP concentrations are interpreted as reference-dose amplitudes. For multi-duration calibration, each scenario keeps the CAP window instantaneous and scales the concentration amplitude by `duration_min / CAP_reference_duration_min`.

- `H2O2/dissipation_coefficient`
- `NO2_/dissipation_coefficient`
- baseline apoptosis, growth, and division settings
- oxygen-threshold settings not explicitly listed above
- domain geometry and boundaries
- `cancer_cell/initial_population`
- `cancer_cell/initial_population/pattern/sphere/radius`

## How `make prepare_calibration_input` transforms the base input

The generator reads `input_inVitro.csv` and writes `input_inVitro_calibration.csv`.

It preserves the same parameter set as the simplified base file, while changing calibration-specific run settings:

- `number_of_time_steps = 7200`
- `time_step = 0.01 h`
- `statistics_interval = 1`
- `diffusion_grid/spatial_resolution = 40`
- `cancer_cell/initial_population = 30000`
- `cancer_cell/initial_population/pattern/sphere/radius = 220.0`
- `output_directory = results_inVitro_calibration`
- `simulation_title = CAP_inVitro_calibration`

It also updates the CAP timing fields consistently:

- `CAP/start_step = 0`
- `CAP/duration_steps = 1`
- `CAP/start_time_h = 0.0`
- `CAP/duration_h = 0.01`

### Auto-scaling of the initial spheroid radius

The generator currently auto-scales:

- `cancer_cell/initial_population/pattern/sphere/radius`

when the calibration template increases the initial cancer-cell count.

The reason is to preserve approximately comparable initial packing density instead of forcing a much larger population into the original small initialization sphere.

For the current workflow, this is the recommended behavior and should stay enabled unless you intentionally want a different initial crowding regime.

## How the experimental Excel data are used

The calibration uses:

`libs/experimental_data/Data Modeling CCA+PDAC 0 to 72h post Gorjet - WITH MEAN 4 experiments.xlsx`

The reduced and ABM workflows both parse the workbook and extract the relevant viability trajectory for the selected cell line and CAP duration.

The comparison logic is:

- choose the requested cell line and nearest available CAP exposure duration
- build the target viability curve over `0 h`, `24 h`, `48 h`, `72 h`
- run the model
- convert live-cell counts to viability relative to the initial model count
- compute a weighted mismatch using the experimental standard deviations

If a trial simulation fails or violates numerical stability, the objective returns a large penalty.

## Command order and what each command does

### 1) `make build`

```bash
make build
```

What it does:

- rebuilds the `ABM4bio` executable
- should be run after code changes
- does not run a simulation by itself

### 2) `make invitro`

```bash
make invitro OMP_THREADS=8
```

What it does:

- runs one BioDynaMo simulation using `input_inVitro.csv`
- writes outputs to `results_inVitro`
- is useful for debugging baseline biology before calibration

### 3) `make prepare_calibration_input`

```bash
make prepare_calibration_input
```

What it does:

- reads `input_inVitro.csv`
- writes `input_inVitro_calibration.csv`
- does not run BioDynaMo
- prepares the longer calibration template with the larger initial population

Timing note:

- CAP is applied as an instantaneous one-step event at `t = 0`
- the configured exposure duration in minutes scales the injected RONS concentration amplitude
- therefore both are written:
  - `CAP/duration_h = 0.01`
  - `CAP/duration_steps = 1`

### 4) `make invitro_calibration`

```bash
make invitro_calibration OMP_THREADS=100

```

What it does:

- runs one BioDynaMo simulation using `input_inVitro_calibration.csv`
- writes outputs to `results_inVitro_calibration`
- does not optimize parameters
- is useful for checking numerical stability and biological plausibility of the calibration template

### 5) `make calibrate_egi1`

```bash
make calibrate_egi1
```

What it does:

- runs the reduced calibration model only
- does not call the BioDynaMo ABM
- keeps the CAP composition fixed at the template values
- compares predicted viability to the experimental `EGI1` trajectory
- is useful as a quick pre-fit sanity check for intracellular response timescales and damage-response parameters

### 6) `make calibrate_egi1_abm`

```bash
make calibrate_egi1_abm OMP_THREADS=100 ABM_DE_MAXITER=1 ABM_DA_MAXITER=10 ABM_NM_MAXITER=40
```

What it does:

- runs ABM-in-the-loop calibration
- repeatedly launches BioDynaMo using temporary copies of `input_inVitro_calibration.csv`
- changes only the selected fit parameters for each trial
- reads each trial `stats.csv`
- converts live cancer-cell counts into simulated viability
- compares simulated viability against the experimental curve
- treats CAP composition as fixed from `input_inVitro_calibration.csv`
- penalizes truncated runs that do not reach the target time horizon
- penalizes late extinction after 24 h / 48 h so the optimizer does not accept collapse to zero just because late experimental SD is large
- writes `calibration_outputs/EGI1_abm/abm_progress.csv` after every evaluation
- uses a lower-evaluation default optimizer sequence (`dual_annealing`, `powell`, then `nelder_mead`) instead of relying on a broad global box search
- if total cells exceed the initial seeded population by even one agent, all positive-time viability predictions for that evaluation are set to `0%`
- selects the parameter set with the lowest weighted mismatch

### Culture scoring modes

The ABM calibration can score the experimental time points in three ways:

- `ABM_CULTURE_SCORE_MODE=continuous`
  - treats 0/24/48/72 h as one continuous simulated culture
  - viability is normalized to the simulated 0 h count

- `ABM_CULTURE_SCORE_MODE=separate`
  - treats terminal endpoint cultures as separately seeded from the same starting population
  - viability is normalized to the intended seeded population

- `ABM_CULTURE_SCORE_MODE=hybrid`
  - averages continuous and separate-culture scoring
  - useful for diagnostic comparison, but not the preferred next calibration route

### What `prepare_calibration_input` does vs calibration make targets

`make prepare_calibration_input` writes **`input_inVitro_calibration.csv` once**. That file is a **template** (default: 72 h horizon, 5 min CAP reference concentrations, cell count, diffusion grid). It is **not** updated when you run `make calibrate_egi1_abm_0to24h_15s` or other calibration targets.

Each calibration make target only changes **command-line flags** passed to `calibrate_cap_viability_abm.py`:

| Make variable | Effect |
|---|---|
| `CAP_DURATIONS_MIN` | Which experimental CAP exposure(s) to fit (e.g. `15/60` = 15 s) |
| `ABM_TARGET_END_TIME_H` | Last time point in the Excel trajectory used for scoring (24, 48, or 72) |
| `ABM_CULTURE_SCORE_MODE` | How viability is scored (`continuous` for 0→endpoint trajectories) |
| `ABM_OUT_DIR` | Where reports and `abm_runs/eval_*/.../input.csv` are written |

For every optimizer trial, the script **copies the template** and writes a fresh `input.csv` under `abm_runs/eval_XXXX/dur_15s/rep_YY/` with:

- `number_of_time_steps` = endpoint hours / `time_step` (e.g. 2400 steps for 24 h at 0.01 h)
- `CAP/duration_h` = physical CAP length in hours (15 s → `15/3600` h)
- fitted parameters from the current candidate

So **`input_inVitro_calibration.csv` stays on disk unchanged**; per-run CAP duration and horizon live in the generated eval inputs.

### Trajectory calibration per CAP duration (preferred)

Use **continuous** scoring over all experimental times from 0 h through the endpoint:

- `make calibrate_egi1_abm_0to24h_15s` — 15 s CAP, simulate and score **0, 24 h** (and intermediate points if present)
- `make calibrate_egi1_abm_0to48h_15s` — 15 s CAP, **0–48 h**
- `make calibrate_egi1_abm_0to72h_15s` — 15 s CAP, **0–72 h**

Repeat for `30s`, `1m`, … `5m` via `calibrate_egi1_abm_0to{24,48,72}h_*`.

Suggested order: calibrate **0–24 h** per duration first, then **0–48 h**, then **0–72 h** (optionally warm-start later stages from the previous stage’s `abm_calibration_report.json`).

Legacy **separate-culture** single-endpoint targets (`calibrate_egi1_abm_24h_15s`, etc.) were removed from the Makefile; use the `0to*` targets above instead.

## Recommended practical order

```bash
source /home/aiwsif/Desktop/ABM4bio/libs/biodynamo-v1.05.143/bin/thisbdm.sh
make show-config
make build
make invitro OMP_THREADS=100

make prepare_calibration_input
make invitro_calibration OMP_THREADS=100 SEED=1234
make calibrate_egi1
make calibrate_egi1_abm_0to24h_15s OMP_THREADS=100 SEED=1234
make calibrate_egi1_abm_0to48h_15s OMP_THREADS=100 SEED=1234
make calibrate_egi1_abm_0to72h_15s OMP_THREADS=100 SEED=1234
```


For more exploratory version use:

```bash
make calibrate_egi1_abm ABM_METHODS="differential_evolution dual_annealing powell nelder_mead"
```

After reviewing the `24 h` best fit, proceed to:

```bash
make clean_calibration
make calibrate_egi1_abm_48h OMP_THREADS=100 SEED=1234

make clean_calibration
make calibrate_egi1_abm_72h OMP_THREADS=100 SEED=1234
```

## Notes and failure handling

- `make prepare_calibration_input` only prepares the template CSV
- `make invitro_calibration` is the fastest one-command test of the calibration template
- `make calibrate_egi1_abm` is the full optimizer-driven ABM loop
- if calibration returns `WRMSE = 1000000`, inspect:
  - `calibration_outputs/EGI1_abm/abm_runs/`
- this usually indicates simulation failure, a truncated run that ended before 72 h, overgrowth beyond the enforced cell limit, or numerical instability during one or more trials
- if most runs terminate near `0.01 h`, the search box is likely entering numerically unstable regions too often; prefer template-seeded local methods and narrower parameter bounds
