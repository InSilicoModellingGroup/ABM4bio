# CAP Calibration with Instantaneous Dose Scaling

> **ABM calibration update:** `calibrate_cap_viability_abm.py` now sets `CAP/duration_h` to the
> physical exposure length in hours (e.g. 15 s → `15/3600` h) and does **not** scale RONS
> concentrations by exposure ratio. Sub-minute exposures still overlap only part of the first
> `time_step` (0.01 h = 36 s); see `describe_cap_abm_timestep()` in the calibration log.

## Problem Statement

The ABM4bio model has a fundamental **timing resolution issue** that prevents realistic simulation of short CAP (cold atmospheric plasma) exposures:

- **Current time step**: `time_step = 0.01 h = 36 seconds`
- **CAP exposures to model**: 15 s, 30 s, 1 min, 2 min, 3 min, 4 min, 5 min

### The Timing Gap

With `dt = 0.01 h`:
| CAP exposure | Equivalent in hours | Required steps | Actual steps |
|---|---|---|---|
| 15 s | 0.00417 h | 0.42 | ⌈0.42⌉ = 1 |
| 30 s | 0.00833 h | 0.83 | ⌈0.83⌉ = 1 |
| 1 min | 0.01667 h | 1.67 | ⌈1.67⌉ = 2 |
| 5 min | 0.08333 h | 8.33 | ⌈8.33⌉ = 9 |

**Result**: 15 s and 30 s exposures **cannot be resolved** with the current timestep.

### The Geometric Problem

Additionally, 30,000 cells in a 25 µm radius sphere is **not geometrically realistic**:
- Cell diameter ≈ 10–20 µm
- Sphere volume ≈ 65,000 µm³
- Space for ~1–5 cells realistically

**Action**: Geometry and baseline proliferation must be fixed independently of CAP calibration.

### The Diffusion-Coefficient Interpretation

The current CAP example keeps:

- `H2O2/diffusion_coefficient = 900`
- `NO2_/diffusion_coefficient = 900`

These RONS diffusion coefficients are treated as effective numerical diffusion parameters, not direct aqueous molecular diffusion coefficients.

Small molecules in water can diffuse on the order of `10^-9 m^2/s`, which is approximately `3.6e6 µm^2/h`. Moving directly to physical diffusion values such as `4e6–7e6 µm^2/h` would make the explicit diffusion stability limit much more restrictive.

For now, keep the effective numerical values. Later, if a modular CAP-specific implicit solver or diffusion substepping path is implemented, the model can move toward physical diffusion coefficients without forcing that numerical burden onto the default workflow.

---

## New Approach: Instantaneous Dose with Concentration Scaling

Instead of trying to discretize sub-minute exposures into ABM timesteps, we treat **CAP as an instantaneous event** at `t = 0` and **scale the initial RONS concentration** by exposure time:

### Key Concepts

1. **CAP is instantaneous**: 
   - `CAP/duration_steps = 1` (always)
   - `CAP/duration_h ≈ time_step` (minimal duration, effectively instantaneous)

2. **RONS dose scales linearly with exposure time**:
   $$C_{\text{H}_2\text{O}_2}(t_{\text{CAP}}) = A_{\text{H}_2\text{O}_2} \cdot t_{\text{CAP}}$$
   $$C_{\text{NO}_2^-}(t_{\text{CAP}}) = A_{\text{NO}_2^-} \cdot t_{\text{CAP}}$$
   where:
   - $t_{\text{CAP}}$ = exposure duration (minutes)
   - $A_{\text{H}_2\text{O}_2, \text{NO}_2^-}$ = amplitude (concentration per minute)

3. **Reference duration**:
   - Default reference: **5 minutes**
   - At 5 min: $C_{\text{H}_2\text{O}_2} = 0.7$, $C_{\text{NO}_2^-} = 0.3$
   - Amplitudes: $A_{\text{H}_2\text{O}_2} = 0.14$/min, $A_{\text{NO}_2^-} = 0.06$/min

### Scaling Examples

Using default reference (5 min) with base concentrations (0.7, 0.3):

| Exposure time | Time ratio | H₂O₂ concentration | NO₂⁻ concentration | Sum |
|---|---|---|---|---|
| 15 s (0.25 min) | 0.05 | 0.035 | 0.015 | 0.050 |
| 30 s (0.5 min) | 0.10 | 0.070 | 0.030 | 0.100 |
| 1 min | 0.20 | 0.140 | 0.060 | 0.200 |
| 2 min | 0.40 | 0.280 | 0.120 | 0.400 |
| 5 min | 1.00 | 0.700 | 0.300 | 1.000 |

---

## Updated Workflow

### 1. Prepare Input for Each CAP Exposure

```bash
# For 15 s exposure (0.25 min)
python3 scripts/prepare_calibration_input.py \
  --cap-exposure-min 0.25 \
  --output input_CAP_0.25min.csv

# For 30 s exposure (0.5 min)
python3 scripts/prepare_calibration_input.py \
  --cap-exposure-min 0.5 \
  --output input_CAP_0.5min.csv

# For 5 min exposure (default reference)
python3 scripts/prepare_calibration_input.py \
  --cap-exposure-min 5 \
  --output input_CAP_5min.csv
```

### 2. Run ABM Simulations

For each input file, run ABM simulations at the corresponding exposure times (15 s, 30 s, 5 min, etc.):

```bash
# Example: 15 s exposure
./build/ABM4bio --config input_CAP_0.25min.csv --seed 1
```

### 3. Compare with Experimental Data

- Extract simulated viability at 24, 48, 72 hours post-exposure
- Compare with experimental measurements
- Adjust biological parameters (_not_ timing or geometry)

### 4. Calibration Parameters (Keep Fixed)

These parameters should **NOT be calibrated** because they fix the structural/timing issues:

| Parameter | Value | Reason |
|---|---|---|
| `time_step` | 0.01 h | Fixes temporal resolution |
| `CAP/duration_steps` | 1 | Instantaneous application |
| `CAP/duration_h` | 0.01 | One timestep |
| `initial_cells` (spheroid) | 30,000 | Geometrically consistent with the updated spheroid radius |
| `sphere/radius` | 220 µm | Keeps cell density realistic |
| `H2O2/diffusion_coefficient` | 900 | Effective numerical diffusion |
| `NO2_/diffusion_coefficient` | 900 | Effective numerical diffusion |

### 5. Calibration Parameters (Can Tune)

These biological parameters control CAP response:

| Parameter | Purpose |
|---|---|
| `CAP/H2O2/concentration` | H2O2 reference-dose amplitude |
| `CAP/NO2_/concentration` | NO2_ reference-dose amplitude |
| `cancer_cell/intracellular/uptake/H2O2` | Intracellular H2O2 uptake |
| `cancer_cell/intracellular/uptake/NO2_` | Intracellular NO2_ uptake |
| `cancer_cell/intracellular/antioxidant/k_scavenge` | Antioxidant scavenging effectiveness |
| `cancer_cell/intracellular/antioxidant/capacity` | Antioxidant buffering capacity |
| `cancer_cell/intracellular/alpha/H2O2` | H2O2 damage weighting |
| `cancer_cell/intracellular/alpha/NO2_` | NO2_ damage weighting |
| `cancer_cell/intracellular/damage/k_induction` | ROS/RNS damage induction rate |
| `cancer_cell/intracellular/damage/threshold` | Damage needed to trigger response |
| `cancer_cell/intracellular/damage/k_repair` | Repair rate |
| `cancer_cell/intracellular/damage/probability` | Probability of damage-mediated response |
| `cancer_cell/can_divide/CAP_sensitivity` | CAP-induced division increase (β_CAP) |

---

## Usage: prepare_calibration_input.py

### New Arguments

```
--cap-exposure-min DURATION
  CAP exposure duration in minutes. 
  Used to scale initial RONS concentrations. 
  CAP applied instantaneously at t=0.
  Default: 5.0 min

--cap-reference-duration-min DURATION
  Reference CAP exposure duration (in minutes) at which 
  the base concentrations are defined (default: 5 min = template reference).
  Default: 5.0 min

--cap-h2o2-base-concentration CONC
  Base CAP H2O2 concentration amplitude at reference duration. 
  Scaled by exposure_min/reference_duration.
  Default: 0.7

--cap-no2-base-concentration CONC
  Base CAP NO2_ concentration amplitude at reference duration. 
  Scaled by exposure_min/reference_duration.
  Default: 0.3
```

### Example Makefile Usage

```makefile
CAP_DURATIONS_MIN ?= 0.25 0.5 1 2 5

prepare_all_cap_inputs:
	@for dur in $(CAP_DURATIONS_MIN); do \
	  python3 scripts/prepare_calibration_input.py \
	    --cap-exposure-min $$dur \
	    --output input_CAP_$${dur}min.csv; \
	done
	@echo "Generated input files for CAP durations: $(CAP_DURATIONS_MIN) min"
```

---

## Why This Works

1. **Avoids timestep resolution problems**: No need to discretize <36 s into timesteps
2. **Linear dose-response**: Biological data suggests RONS damage ∝ exposure time
3. **Separates concerns**: 
   - Geometry/timing/baseline proliferation → structural fixes
   - CAP-response parameters → calibration targets
4. **Computational efficiency**: 5× faster than simulating duration_steps > 1
5. **Clean parameter space**: Concentration rather than duration is the true biological dose

---

## Calibration Strategy

### Phase 1: Baseline Proliferation Reference

- Establish baseline proliferation with:
  - Correct spheroid size and cell density
  - Correct cell cycle parameters
  - The 15 s CAP condition as the minimum-dose reference

Baseline proliferation was adjusted using the 15 s CAP condition as a minimum-dose reference, not as an untreated biological control.

### Phase 2: Reference Calibration (5 min CAP)

- Fit CAP-response parameters against 5 min exposure data
- Example parameters:
  ```
  cancer_cell/intracellular/damage/threshold = 0.5–2.0
  cancer_cell/intracellular/damage/k_repair = 0.001–0.02
  cancer_cell/intracellular/antioxidant/k_scavenge = 0.01–0.1
  ```

### Phase 3: Duration Scaling Validation

- Apply fitted parameters to 15 s, 30 s, 1 min, 2 min exposures
- Compare predicted viability at 24, 48, 72 h against experimental data
- If dose-response is non-linear, refine damage accumulation model

---

## Notes on Nomenclature

The **15 s CAP exposure** is the **minimum-dose reference**, not an untreated control:
- 15 s still applies CAP → cells still receive RONS damage
- Use **0 min (no exposure)** as the true control for proliferation baseline
- Report as: "minimum-dose reference" to avoid confusion

---

## Troubleshooting

### Problem: Scaled concentrations exceed 1.0

**Cause**: Exposure time × base amplitude > 1.0 (physically unrealistic)

**Solution**:
1. Reduce base concentrations: `--cap-h2o2-base-concentration 0.35 --cap-no2-base-concentration 0.15`
2. Change reference duration: `--cap-reference-duration-min 10` (makes longer exposures "normal")
3. Check experimental setup (maybe very high RONS delivery at 5 min?)

### Problem: Model doesn't show dose-response

**Causes**:
- Damage threshold too high → cells ignore RONS
- Repair rate too fast → cells heal faster than damage accumulates
- Concentration scaling factor wrong

**Debug steps**:
1. Manually set very high concentrations: `--cap-h2o2-base-concentration 1.0`
2. Confirm model shows increased cell death
3. Then calibrate threshold and repair rate
