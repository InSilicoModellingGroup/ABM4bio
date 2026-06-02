#!/usr/bin/env python3
"""Report viable cancer cells at standard time points from stats.csv."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from abm_io import viable_counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stats_csv", type=Path, nargs="?", default=Path("results_control/stats.csv"))
    parser.add_argument(
        "--times-h",
        type=float,
        nargs="+",
        default=[0.0, 1.0, 2.0, 4.0, 6.0, 12.0, 24.0, 48.0, 72.0],
    )
    args = parser.parse_args()

    stats = pd.read_csv(args.stats_csv, skipinitialspace=True)
    times = stats["current_time"].to_numpy(dtype=float)
    viable = viable_counts(stats)
    total = (
        stats["N_cells"].to_numpy(dtype=float)
        if "N_cells" in stats.columns
        else viable
    )

    print(f"stats: {args.stats_csv.resolve()}")
    print(f"frames: {len(stats)}  t_end: {times[-1]:.2f} h")
    print("time_h  viable  N_total  pct_vs_t0")
    baseline = float(viable[0])
    for target_h in args.times_h:
        idx = int(np.argmin(np.abs(times - target_h)))
        v = float(viable[idx])
        n = float(total[idx])
        t = float(times[idx])
        pct = 100.0 * v / max(baseline, 1.0e-12)
        print(f"{t:6.2f}  {v:8.1f}  {n:8.1f}  {pct:7.1f}%")


if __name__ == "__main__":
    main()
