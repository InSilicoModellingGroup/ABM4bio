# Oxygen Gradient Diffusion Model

## Overview
This example demonstrates a simple oxygen gradient diffusion model where cells respond to oxygen concentration gradients through chemotaxis, growth, division, and apoptosis.

## Model Description

### Key Features
- **Single diffusive substance**: Oxygen (O2) is the only biochemical modeled
- **Positive chemotaxis**: Cells migrate towards higher oxygen concentrations
- **Oxygen-dependent behaviors**:
  - **Growth**: Cells grow faster in higher oxygen environments
  - **Division**: Cells require sufficient oxygen to divide
  - **Apoptosis**: Cells undergo apoptosis in low oxygen conditions

### Parameters
- **Simulation domain**: 3D bounded polar domain (-100 to +100 units)
- **Initial cell population**: 100 cells in a spherical cluster (radius 5.0)
- **Oxygen gradient**: Initial values range from 1.0 to 5.0
- **Oxygen diffusion coefficient**: 1.0e+2
- **Oxygen dissipation**: Small dissipation (0.01) to create gradients

### Cell Behaviors

#### Migration (Chemotaxis)
- Cells migrate towards **higher oxygen concentrations** (positive chemotaxis coefficient: 0.001)
- High migration probability (0.5) with gradient normalization enabled
- Oxygen threshold for chemotaxis: 0.5

#### Growth
- Cells grow when oxygen levels are above 1.0
- Oxygen-dependent growth rate: 1.0

#### Division
- Cells divide when they reach diameter cutoff (14.5)
- Requires oxygen levels above 2.0
- Maximum 150 divisions per cell

#### Apoptosis
- Cells undergo apoptosis when oxygen drops below 0.5
- High apoptosis probability (0.8) in low oxygen conditions

#### Oxygen Consumption
- Cells consume oxygen at a net rate of -0.01
- Small variation in consumption (std: 0.001)

#### Optional: Michaelis-Menten Kinetics
By default, substance consumption/production uses a constant rate (`net_balance`).
You can switch to **Michaelis-Menten kinetics** for any substance by adding these
CSV parameters (shown here for O2, but applicable to any substance):

```
cancer_cell/O2/secretion/kinetics_model,string,michaelis_menten
cancer_cell/O2/secretion/michaelis_menten/Vmax,float,0.005
cancer_cell/O2/secretion/michaelis_menten/Km,float,0.5
```

| Parameter | Description |
|-----------|-------------|
| `kinetics_model` | `constant` (default) or `michaelis_menten` |
| `michaelis_menten/Vmax` | Maximum consumption/production rate (must be > 0) |
| `michaelis_menten/Km` | Half-saturation constant — the concentration at which the rate is 50% of Vmax (must be > 0) |

The effective rate becomes: **R = ±Vmax × C / (Km + C)**, where:
- The sign follows `net_balance` (negative = consumption, positive = production)
- `C` is the local substrate concentration
- The existing `net_balance/std` stochastic variation is also applied
- Saturation and dependency settings still function normally

When `kinetics_model` is omitted or set to `constant`, the original flat-rate
model is used unchanged (full backward compatibility).

## Running the Simulation

### Build the project
```bash
make build
```

### Run the simulation
```bash
make invitro
```

### Analyze results
```bash
make analyze
```

### Clean up
```bash
make clean
```

## Expected Results
- Cells initially clustered in the center will migrate towards regions with higher oxygen concentration
- Cell population dynamics will reflect the balance between:
  - Growth and division in oxygen-rich regions
  - Apoptosis in oxygen-depleted regions
  - Oxygen consumption by cells

## Output Files
- `results_inVitro/`: Main simulation output directory
  - Paraview visualization files
  - `stats.csv`: Time series data
- `oxygen_gradient_analysis.png`: Plots of cell population dynamics
- `population_summary.txt`: Summary statistics

## Visualization
Use Paraview to visualize the 3D cell distribution and oxygen gradient over time.
