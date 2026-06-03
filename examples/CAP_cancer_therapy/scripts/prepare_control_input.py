#!/usr/bin/env python3
"""Build a fast untreated-control ABM input (0–72 h, no CAP, coarse diffusion grid)."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from abm_io import get_float_param, get_param, read_rows, set_param, write_rows  # noqa: F401 (get_float_param used in print)


def scale_sphere_radius(rows: list[list[str]], old_cells: int, new_cells: int) -> None:
    if old_cells <= 0 or new_cells <= 0:
        return
    key = "cancer_cell/initial_population/pattern/sphere/radius"
    raw = get_param(rows, key)
    if raw is None:
        return
    radius = float(raw) * (float(new_cells) / float(old_cells)) ** (1.0 / 3.0)
    set_param(rows, key, f"{radius:.12g}")


def disable_cap(rows: list[list[str]]) -> None:
    set_param(rows, "CAP/enabled", "false", "bool")
    set_param(rows, "CAP/duration_h", "0")
    set_param(rows, "CAP/duration_steps", "0", "int")
    set_param(rows, "CAP/H2O2/concentration", "0")
    set_param(rows, "CAP/NO2_/concentration", "0")


def neutralize_ros(rows: list[list[str]]) -> None:
    """Keep biochemical grids for ABM4bio, but remove ROS coupling for untreated growth."""
    for key, value in (
        ("cancer_cell/intracellular/uptake/H2O2", "0"),
        ("cancer_cell/intracellular/uptake/NO2_", "0"),
        ("cancer_cell/intracellular/damage/k_induction", "0"),
        ("cancer_cell/intracellular/damage/probability", "0"),
        ("cancer_cell/can_divide/CAP_sensitivity", "0"),
    ):
        set_param(rows, key, value)
    for species in ("H2O2", "NO2_"):
        set_param(rows, f"{species}/initial_value", "0", "float")
        set_param(rows, f"{species}/diffusion_coefficient", "0", "float")
        set_param(rows, f"{species}/dissipation_coefficient", "0", "float")


def apply_growth_defaults(rows: list[list[str]]) -> None:
    """ABM4bio requires full mechanism blocks when can_apoptose / can_divide are enabled.

    Defaults target net proliferation (Excel EGI1 control ~2x @ 24 h). Previous apoptose=0.08
    collapsed the culture within hours (8% roll per mechanism check).
    """
    set_param(rows, "cancer_cell/can_apoptose", "true", "bool")
    set_param(rows, "cancer_cell/can_apoptose/probability", "0.0029217516", "float")
    set_param(rows, "cancer_cell/can_apoptose/probability_increment_with_age", "0.0", "float")
    set_param(rows, "cancer_cell/can_apoptose/time_window", "2500", "int")
    set_param(rows, "cancer_cell/can_apoptose/time_window/to_delete", "5000", "int")
    set_param(rows, "cancer_cell/can_grow/probability", "0.52025049", "float")
    set_param(rows, "cancer_cell/can_grow/diameter_rate", "0.4556126", "float")
    set_param(rows, "cancer_cell/can_divide/probability", "0.84592981", "float")
    set_param(rows, "cancer_cell/can_divide/probability_increment_with_age", "0.0", "float")
    set_param(rows, "cancer_cell/can_divide/time_window", "319", "int")
    set_param(rows, "cancer_cell/principal/min", "1", "float")
    set_param(rows, "cancer_cell/principal/max", "1", "float")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-input", type=Path, default=Path("input_control_base.csv"))
    parser.add_argument("--output", type=Path, default=Path("input_control.csv"))
    parser.add_argument("--results-directory", default="results_control")
    parser.add_argument("--simulation-title", default="control_0to72h")
    parser.add_argument("--initial-cells", type=int, default=100)
    parser.add_argument("--total-cell-limit", type=int, default=40000)
    parser.add_argument("--total-hours", type=float, default=72.0)
    parser.add_argument("--time-step-h", type=float, default=0.01)
    parser.add_argument(
        "--grid-resolution",
        type=int,
        default=10,
        help="Diffusion grid points per axis (ABM4bio minimum 10, must be even).",
    )
    parser.add_argument(
        "--viz-interval-h",
        type=float,
        default=1.0,
        help="Hours between ParaView VTK outputs (default: every 1 h).",
    )
    args = parser.parse_args()

    rows = read_rows(args.base_input)
    old_cells = int(round(get_float_param(rows, "cancer_cell/initial_population", args.initial_cells)))
    n_steps = int(round(args.total_hours / args.time_step_h))
    viz_steps = max(1, int(round(args.viz_interval_h / args.time_step_h)))

    set_param(rows, "output_directory", args.results_directory)
    set_param(rows, "simulation_title", args.simulation_title)
    set_param(rows, "number_of_time_steps", str(n_steps), "int")
    set_param(rows, "time_step", f"{args.time_step_h:g}")
    set_param(rows, "statistics_interval", "1", "int")
    set_param(rows, "visualization_interval", str(viz_steps), "int")
    set_param(rows, "simulation/early_stop_on_total_cells_exceeded", "true", "bool")
    set_param(rows, "simulation/total_cell_limit", str(int(args.total_cell_limit)), "int")
    grid_res = int(args.grid_resolution)
    if grid_res < 10:
        grid_res = 10
    if grid_res % 2 != 0:
        grid_res -= 1
    set_param(rows, "diffusion_grid/spatial_resolution", str(grid_res), "int")
    set_param(rows, "cancer_cell/initial_population", str(args.initial_cells), "int")
    scale_sphere_radius(rows, old_cells, args.initial_cells)

    disable_cap(rows)
    neutralize_ros(rows)
    apply_growth_defaults(rows)

    write_rows(args.output, rows)
    print(f"Wrote {args.output}")
    print(f"  steps: {n_steps} ({args.total_hours:g} h @ {args.time_step_h:g} h/step)")
    print(f"  grid resolution: {get_param(rows, 'diffusion_grid/spatial_resolution')}")
    print(f"  initial cells: {args.initial_cells}")
    print(f"  total cell limit: {int(args.total_cell_limit)}")
    print(f"  visualization_interval: {viz_steps} steps (~{args.viz_interval_h:g} h, ~{n_steps // viz_steps + 1} frames)")
    print(
        f"  growth (LM template): "
        f"apoptose_p={get_float_param(rows, 'cancer_cell/can_apoptose/probability'):.7g}  "
        f"divide_p={get_float_param(rows, 'cancer_cell/can_divide/probability'):.7g}  "
        f"grow_rate={get_float_param(rows, 'cancer_cell/can_grow/diameter_rate'):.7g}  "
        f"divide_window={int(get_float_param(rows, 'cancer_cell/can_divide/time_window'))}"
    )
    print("  CAP: disabled | ROS uptake/damage: off")


if __name__ == "__main__":
    main()
