"""ABM input CSV helpers, stats parsing, and viability scoring."""

from __future__ import annotations

import csv
import math
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROBABILITY_BOUNDS = (0.01, 0.99)  # hard cap: probabilities must stay in (0, 1)
PROBABILITY_MAX = PROBABILITY_BOUNDS[1]  # 0.99 — never allow >= 1.0


@dataclass(frozen=True)
class FitParameterization:
    names: Tuple[str, ...]
    bounds: Tuple[Tuple[float, float], ...]


# Default 4-parameter set derived from LM calibration baseline (input_control.csv):
#   can_divide/probability   LM value: 0.8459  bounds: [0.30, PROBABILITY_MAX]
#   can_apoptose/probability LM value: 0.0029  bounds: [0.0001, 0.10]  (low but non-trivial)
#   can_grow/diameter_rate   LM value: 0.4556  bounds: [0.05, 1.5]     (timescale-sensitive, NOT a probability)
#   can_grow/probability     LM value: 0.5203  bounds: [0.01, PROBABILITY_MAX]
# Fixed at LM template (not calibrated here):
#   can_divide/time_window   LM value: 319     (structural, set by LM stage)
CONTROL_GROWTH_PARAMS = FitParameterization(
    names=(
        "cancer_cell/can_divide/probability",
        "cancer_cell/can_apoptose/probability",
        "cancer_cell/can_grow/diameter_rate",
        "cancer_cell/can_grow/probability",
        "cancer_cell/can_divide/time_window",
    ),
    bounds=(
        (0.30, PROBABILITY_MAX),       # can_divide/probability: never >= 1
        (0.0001, 0.25),                # can_apoptose/probability: never >= 1
        (0.05, 1.5),                   # can_grow/diameter_rate: growth rate, not a probability
        (0.01, 0.99),                  # can_grow/probability: LM=0.52, cap at 0.70 (>0.70 causes extreme overgrowth)
        (50.0, 500.0),                 # can_divide/time_window: integer steps, not a probability
    ),
)


def read_rows(path: Path) -> List[List[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [row for row in csv.reader(handle)]


def write_rows(path: Path, rows: Sequence[Sequence[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerows(rows)


def set_param(rows: List[List[str]], key: str, value: str, param_type: str | None = None) -> None:
    for row in rows:
        if len(row) >= 3 and row[0].strip() == key:
            row[2] = value
            if param_type is not None and len(row) >= 2:
                row[1] = param_type
            return
    rows.append([key, param_type or "string", value])


def get_param(rows: Sequence[Sequence[str]], key: str) -> str | None:
    for row in rows:
        if len(row) >= 3 and row[0].strip() == key:
            return row[2]
    return None


def remove_param(rows: List[List[str]], key: str) -> None:
    rows[:] = [row for row in rows if len(row) < 3 or row[0].strip() != key]


def get_float_param(rows: Sequence[Sequence[str]], key: str, default: float = 0.0) -> float:
    value = get_param(list(rows), key)
    return float(value) if value is not None else default


def select_fit_parameters(selected: Sequence[str]) -> FitParameterization:
    names = tuple(selected)
    if not names:
        raise ValueError("At least one fit parameter is required.")
    lookup = {name: index for index, name in enumerate(CONTROL_GROWTH_PARAMS.names)}
    unknown = [name for name in names if name not in lookup]
    if unknown:
        raise ValueError(f"Unknown fit parameters: {unknown}")
    return FitParameterization(
        names=names,
        bounds=tuple(CONTROL_GROWTH_PARAMS.bounds[lookup[name]] for name in names),
    )


def clip_to_bounds(x: np.ndarray, bounds: Sequence[Tuple[float, float]]) -> np.ndarray:
    clipped = np.asarray(x, dtype=float).copy()
    for index, (lo, hi) in enumerate(bounds):
        clipped[index] = min(hi, max(lo, clipped[index]))
    return clipped


def template_point(rows: Sequence[Sequence[str]], fit: FitParameterization) -> np.ndarray:
    return np.asarray([get_float_param(rows, name, 0.0) for name in fit.names], dtype=float)


def viable_counts(stats: pd.DataFrame, phenotype_id: int = 1) -> np.ndarray:
    stats = stats.copy()
    stats.columns = [column.strip() for column in stats.columns]
    base = f"N_cells_pheno_{phenotype_id}"
    phase_cols = [f"{base}_{phase}" for phase in ("G1", "Sy", "G2", "Di", "Tr")]
    if all(column in stats.columns for column in phase_cols):
        total = np.zeros(len(stats), dtype=float)
        for column in phase_cols:
            total += stats[column].to_numpy(dtype=float)
        return np.maximum(0.0, total)
    apoptotic_col = f"{base}_Ap"
    if base in stats.columns and apoptotic_col in stats.columns:
        return np.maximum(
            0.0,
            stats[base].to_numpy(dtype=float) - stats[apoptotic_col].to_numpy(dtype=float),
        )
    if "N_cells" in stats.columns:
        return stats["N_cells"].to_numpy(dtype=float)
    raise ValueError("stats.csv does not contain recognizable cancer-cell count columns.")


def population_bounds(
    target_curve: pd.DataFrame,
    initial_population: int,
    std_scale: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    curve = target_curve.sort_values("time_h").reset_index(drop=True)
    times = curve["time_h"].to_numpy(dtype=float)
    mean_signal = curve["mean_signal"].to_numpy(dtype=float)
    sd_signal = (
        curve["sd_signal"].to_numpy(dtype=float)
        if "sd_signal" in curve.columns
        else np.zeros_like(mean_signal)
    )
    sd_signal = np.where(np.isfinite(sd_signal), sd_signal, 0.0)
    baseline_idx = int(np.argmin(np.abs(times)))
    baseline = float(mean_signal[baseline_idx])
    if not np.isfinite(baseline) or baseline <= 0.0:
        default = np.asarray([float(max(initial_population, 1))], dtype=float)
        return np.asarray([0.0]), np.asarray([0.0]), default

    upper = np.maximum(0.0, mean_signal + std_scale * sd_signal)
    lower = np.maximum(0.0, mean_signal - std_scale * sd_signal)
    ratio_upper = np.maximum(upper / baseline, 0.0)
    ratio_lower = np.maximum(lower / baseline, 0.0)
    population0 = float(max(initial_population, 1))
    ceiling = np.maximum(population0 * ratio_upper, population0 * ratio_lower)
    floor = np.maximum(population0 * ratio_lower, 0.0)
    ceiling = np.maximum(ceiling, 1.0)
    return times.astype(float), floor.astype(float), ceiling.astype(float)


def evaluate_growth_trace(
    times: np.ndarray,
    viable: np.ndarray,
    total: np.ndarray,
    target_curve: pd.DataFrame,
    target_times: np.ndarray,
    required_end_h: float,
    first_alive_h: float,
    initial_population: int,
    floor_times: np.ndarray,
    floor_counts: np.ndarray,
    ceiling_times: np.ndarray,
    ceiling_counts: np.ndarray,
    *,
    truncation_scale: float = 12.0,
    early_extinction_scale: float = 8.0,
    late_extinction_scale: float = 8.0,
    late_viability_floor_pct: float = 25.0,
) -> Tuple[float, np.ndarray, Dict[str, object]]:
    if len(times) == 0 or len(viable) == 0:
        nan = np.full(len(target_times), np.nan)
        return 1.0e6, nan, {"invalid_metrics": True, "truncated": True}

    predicted_counts = np.interp(target_times, times, viable, left=viable[0], right=viable[-1])
    viability_pct = 100.0 * predicted_counts / max(float(initial_population), 1.0)

    target_viability = target_curve["target_viability_pct"].to_numpy(dtype=float)
    target_sd = target_curve["target_sd_pct"].to_numpy(dtype=float)
    target_sd = np.where(np.isfinite(target_sd) & (target_sd > 1.0e-9), target_sd, 1.0)

    eps = 1.0
    log_pred = np.log(np.maximum(viability_pct[1:], eps))
    log_target = np.log(np.maximum(target_viability[1:], eps))
    log_sd = np.log1p(target_sd[1:] / np.maximum(target_viability[1:], eps))
    log_sd = np.where(log_sd > 1.0e-9, log_sd, 1.0)
    wrmse = float(np.sqrt(np.mean(np.square((log_pred - log_target) / log_sd))))

    final_h = float(times[-1])
    completion = min(1.0, final_h / max(required_end_h, 1.0e-12))
    truncation_penalty = truncation_scale * max(0.0, 1.0 - completion) ** 2

    zero_idx = np.flatnonzero(viable <= 0)
    first_zero_h = float(times[int(zero_idx[0])]) if len(zero_idx) else None
    early_penalty = 0.0
    if first_zero_h is not None and first_zero_h < first_alive_h:
        early_penalty = early_extinction_scale * (
            (first_alive_h - first_zero_h) / max(first_alive_h, 1.0e-12)
        ) ** 2

    late_mask = (target_times >= first_alive_h) & (
        target_viability >= late_viability_floor_pct
    )
    late_shortfall = np.maximum(0.0, late_viability_floor_pct - viability_pct[late_mask])
    late_penalty = (
        late_extinction_scale
        * float(np.mean(np.square(late_shortfall / max(late_viability_floor_pct, 1.0e-12))))
        if np.any(late_mask)
        else 0.0
    )

    sim_floor = np.interp(times, floor_times, floor_counts, left=floor_counts[0], right=floor_counts[-1])
    sim_ceiling = np.interp(
        times, ceiling_times, ceiling_counts, left=ceiling_counts[0], right=ceiling_counts[-1]
    )
    undergrowth = bool(np.any(viable < sim_floor))
    overgrowth = bool(np.any(total > sim_ceiling))
    undergrowth_penalty = 5.0 if undergrowth else 0.0
    if overgrowth:
        max_excess = float(np.max(total / np.maximum(sim_ceiling, 1.0)))
        overgrowth_penalty = 5.0 * max_excess  # proportional: 2x ceiling -> 10, 10x -> 50
    else:
        overgrowth_penalty = 0.0

    score = wrmse + truncation_penalty + early_penalty + late_penalty + undergrowth_penalty + overgrowth_penalty
    if not np.isfinite(score):
        score = 1.0e6

    details = {
        "wrmse": wrmse,
        "truncation_penalty": truncation_penalty,
        "early_extinction_penalty": early_penalty,
        "late_extinction_penalty": late_penalty,
        "undergrowth_penalty": undergrowth_penalty,
        "overgrowth_penalty": overgrowth_penalty,
        "final_time_h": final_h,
        "completion_ratio": completion,
        "first_zero_time_h": first_zero_h,
        "max_total": float(np.max(total)) if len(total) else 0.0,
        "undergrowth_below_floor": undergrowth,
        "overgrowth_above_ceiling": overgrowth,
        "truncated": final_h + 1.0e-9 < required_end_h,
        "invalid_metrics": not np.isfinite(score),
    }
    return score, viability_pct, details


def compact_output_dir(output_dir: Path) -> None:
    if not output_dir.exists():
        return
    for child in output_dir.iterdir():
        if child.name == "stats.csv":
            continue
        if child.is_dir():
            shutil.rmtree(child, ignore_errors=True)
        else:
            try:
                child.unlink()
            except OSError:
                pass


def probe_stats_cells(stats_path: Path) -> Tuple[float, int]:
    """Return (current_time_h, N_cells) from the last non-header line of stats.csv.

    Returns (-1.0, -1) if the file is absent, unreadable, or has no data rows yet.
    """
    if not stats_path.exists() or stats_path.stat().st_size < 20:
        return -1.0, -1
    try:
        with stats_path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            chunk = min(handle.tell(), 8192)
            handle.seek(-chunk, os.SEEK_END)
            tail = handle.read().decode("utf-8", errors="replace")
        # Parse header to locate the N_cells column index
        all_lines = tail.splitlines()
        header_line: Optional[str] = None
        for line in all_lines:
            stripped = line.strip()
            if stripped.startswith("current_time"):
                header_line = stripped
                break
        if header_line is None:
            return -1.0, -1
        headers = [h.strip() for h in header_line.split(",")]
        try:
            n_cells_idx = headers.index("N_cells")
        except ValueError:
            return -1.0, -1
        data_lines = [
            line.strip()
            for line in all_lines
            if line.strip() and not line.strip().startswith("current_time")
        ]
        if not data_lines:
            return -1.0, -1
        parts = data_lines[-1].split(",")
        current_h = float(parts[0])
        n_cells = int(float(parts[n_cells_idx]))
        return current_h, n_cells
    except (OSError, ValueError, IndexError):
        return -1.0, -1


def probe_stats_progress(stats_path: Path, required_end_h: float) -> Tuple[float, float]:
    if not stats_path.exists() or stats_path.stat().st_size < 20:
        return 0.0, 0.0
    try:
        with stats_path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            chunk = min(handle.tell(), 8192)
            handle.seek(-chunk, os.SEEK_END)
            tail = handle.read().decode("utf-8", errors="replace")
        lines = [
            line.strip()
            for line in tail.splitlines()
            if line.strip() and not line.startswith("current_time")
        ]
        if not lines:
            return 0.0, 0.0
        current_h = float(lines[-1].split(",", 1)[0])
        fraction = min(1.0, max(0.0, current_h / max(required_end_h, 1.0e-12)))
        return current_h, fraction
    except (OSError, ValueError, IndexError):
        return 0.0, 0.0


def sensitivity_report(
    objective: object,
    best_x: np.ndarray,
    best_score: float,
    *,
    perturbation: float = 0.05,
) -> List[Dict[str, object]]:
    """Perturb each parameter ±perturbation fraction from best_x and report score delta.

    Also computes a finite-difference Jacobian of the *predictions* (not the scalar score)
    at best_x, then derives:
      - local_sensitivity_a: local curvature proxy from J^T W J (NOT a real parameter
        uncertainty / standard error — for true uncertainty use bootstrap or ABC)
      - corr_a: parameter correlation matrix derived from J^T J (local approximation only)

    Returns a list of dicts (one per parameter) with keys:
        name, best_value, plus_value, minus_value,
        score_plus, score_minus, delta_plus, delta_minus, sensitivity,
        local_sensitivity_a, local_sensitivity_a_pct
    Plus two extra keys in rows[0]: 'correlation_matrix', 'correlation_names'.
    """
    saved_max_evals = objective.max_evals
    objective.max_evals = 0
    objective.set_method("sensitivity")

    n_params = len(objective.fit.names)
    n_times = len(objective.target_times)
    rows: List[Dict[str, object]] = []

    # ---- score-based sensitivity (±perturbation of each param) ----
    J_pred = np.zeros((n_times, n_params))  # Jacobian: dpred/dparam
    pred_plus_all: List[np.ndarray] = []
    pred_minus_all: List[np.ndarray] = []

    for i, name in enumerate(objective.fit.names):
        lo, hi = objective.fit.bounds[i]
        base = float(best_x[i])
        plus_val = min(hi, base * (1.0 + perturbation))
        minus_val = max(lo, base * (1.0 - perturbation))
        step = plus_val - minus_val  # actual step used (central diff)

        x_plus = best_x.copy()
        x_plus[i] = plus_val
        score_plus, pred_plus, _ = objective.simulate(x_plus)

        x_minus = best_x.copy()
        x_minus[i] = minus_val
        score_minus, pred_minus, _ = objective.simulate(x_minus)

        pred_plus_all.append(pred_plus)
        pred_minus_all.append(pred_minus)

        if step > 1e-12 and not np.any(np.isnan(pred_plus)) and not np.any(np.isnan(pred_minus)):
            J_pred[:, i] = (pred_plus - pred_minus) / step

        delta_plus = score_plus - best_score
        delta_minus = score_minus - best_score
        rows.append({
            "name": name,
            "best_value": base,
            "plus_value": plus_val,
            "minus_value": minus_val,
            "score_plus": score_plus,
            "score_minus": score_minus,
            "delta_plus": delta_plus,
            "delta_minus": delta_minus,
            "sensitivity": max(abs(delta_plus), abs(delta_minus)),
            "local_sensitivity_a": float("nan"),
            "local_sensitivity_a_pct": float("nan"),
        })

    # ---- Jacobian-based covariance: J^T W J where W = identity (log-space equal weighting) ----
    try:
        JtJ = J_pred.T @ J_pred  # (n_params x n_params)
        # regularise slightly for near-singular cases
        reg = 1e-8 * np.trace(JtJ) / max(n_params, 1)
        JtJ_reg = JtJ + reg * np.eye(n_params)
        covar = np.linalg.inv(JtJ_reg)
        sigma_a = np.sqrt(np.maximum(np.diag(covar), 0.0))

        # correlation matrix
        outer = np.outer(sigma_a, sigma_a)
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = np.where(outer > 0, covar / outer, 0.0)
        np.fill_diagonal(corr, 1.0)

        for i, row in enumerate(rows):
            row["local_sensitivity_a"] = float(sigma_a[i])
            base = abs(float(row["best_value"]))
            row["local_sensitivity_a_pct"] = float(100.0 * sigma_a[i] / base) if base > 1e-12 else float("nan")

        corr_list = corr.tolist()
        corr_names = [r["name"] for r in rows]
        rows[0]["correlation_matrix"] = corr_list
        rows[0]["correlation_names"] = corr_names
    except (np.linalg.LinAlgError, ValueError):
        pass

    rows.sort(key=lambda r: r["sensitivity"], reverse=True)
    objective.max_evals = saved_max_evals
    return rows


def print_sensitivity_report(rows: List[Dict[str, object]], perturbation: float = 0.05) -> None:
    pct = int(round(perturbation * 100))
    print(f"\n=== Local Parameter Sensitivity (±{pct}% perturbation from best) ===")
    print(f"  NOTE: 'loc_sens' is a local one-at-a-time curvature proxy, NOT a parameter uncertainty estimate.")
    print(f"  {'Parameter':<45} {'Best':>10} {'Δ(+%d%%)' % pct:>10} {'Δ(-%d%%)' % pct:>10}  Sensitivity  loc_sens   ±%")
    print("  " + "-" * 104)
    for r in rows:
        sigma_str = f"{r['local_sensitivity_a']:.4g}" if np.isfinite(float(r["local_sensitivity_a"])) else "  n/a  "
        pct_str = f"{r['local_sensitivity_a_pct']:.1f}%" if np.isfinite(float(r["local_sensitivity_a_pct"])) else "  n/a"
        print(
            f"  {r['name']:<45} {r['best_value']:>10.4g}"
            f"  {r['delta_plus']:>+9.4f}  {r['delta_minus']:>+9.4f}"
            f"  {r['sensitivity']:>10.4f}  {sigma_str:>8}  {pct_str:>6}"
        )
    print()

    # Print correlation matrix if present
    corr = None
    names = None
    for r in rows:
        if "correlation_matrix" in r:
            corr = np.array(r["correlation_matrix"])
            names = r["correlation_names"]
            break
    if corr is not None and names is not None:
        short = [n.split("/")[-1][:12] for n in names]
        print("=== Parameter Correlation Matrix ===")
        header = "  " + " " * 14 + "".join(f"{s:>10}" for s in short)
        print(header)
        for i, row_name in enumerate(short):
            vals = "".join(
                f"  {'1.000':>7}" if j == i else f"  {corr[i, j]:>+7.3f}"
                for j in range(len(short))
            )
            print(f"  {row_name:<14}{vals}")
        print()
        high_pairs = [
            (abs(corr[i, j]), names[i], names[j], corr[i, j])
            for i in range(len(names)) for j in range(i + 1, len(names))
            if abs(corr[i, j]) > 0.8
        ]
        if high_pairs:
            print("  ⚠ Highly correlated pairs (|r| > 0.8):")
            for _, na, nb, val in sorted(high_pairs, reverse=True):
                print(f"    {na.split('/')[-1]} ↔ {nb.split('/')[-1]}  r={val:+.3f}")
            print()


def prediction_at_time(target_times: np.ndarray, preds: np.ndarray, time_h: float) -> float:
    matches = np.flatnonzero(np.isclose(target_times, time_h, atol=1.0e-9))
    if not len(matches):
        return float("nan")
    value = float(preds[int(matches[0])])
    return value if np.isfinite(value) else float("nan")
