#!/usr/bin/env python3
"""Calibrate untreated in-vitro growth (0–72 h) against Excel control data."""

from __future__ import annotations

import argparse
import collections
import json
import math
import os
import signal
import shlex
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Deque, Dict, List, Optional, Sequence, Tuple

# Keep Python-side linear algebra single-threaded by default.
# The ABM subprocess receives its own `OMP_NUM_THREADS` setting separately.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("BLIS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy.optimize import dual_annealing, minimize

try:
    from threadpoolctl import threadpool_limits
except ImportError:
    threadpool_limits = None

# Bayesian Optimization for expensive, noisy ABM simulations
def _ensure_skopt():
    """Try to import skopt, adding common site-packages paths if needed."""
    import sys
    from pathlib import Path
    
    try:
        from skopt import gp_minimize
        from skopt.space import Real
        from skopt.utils import use_named_args
        return True, gp_minimize, Real, use_named_args
    except ImportError:
        # Try adding miniconda site-packages
        candidate_paths = [
            Path.home() / "miniconda" / "lib" / "python3.13" / "site-packages",
            Path.home() / "miniconda3" / "lib" / "python3.13" / "site-packages",
            Path("/usr/local/lib/python3.13/site-packages"),
        ]
        for path in candidate_paths:
            if path.exists() and str(path) not in sys.path:
                sys.path.insert(0, str(path))
        try:
            from skopt import gp_minimize
            from skopt.space import Real
            from skopt.utils import use_named_args
            return True, gp_minimize, Real, use_named_args
        except ImportError:
            return False, None, None, None

_SKOPT_RESULT = _ensure_skopt()
SKOPT_AVAILABLE = _SKOPT_RESULT[0]
if SKOPT_AVAILABLE:
    _gp_minimize = _SKOPT_RESULT[1]
    _Real = _SKOPT_RESULT[2]
    _use_named_args = _SKOPT_RESULT[3]

from abm_io import (
    CONTROL_GROWTH_PARAMS,
    FitParameterization,
    clip_to_bounds,
    compact_output_dir,
    evaluate_growth_trace,
    get_float_param,
    population_bounds,
    prediction_at_time,
    print_sensitivity_report,
    probe_stats_cells,
    probe_stats_progress,
    read_rows,
    select_fit_parameters,
    sensitivity_report,
    set_param,
    template_point,
    viable_counts,
    write_rows,
)
from excel_data import build_control_target, parse_workbook

_PROGRESS_WIDTH = 50
_STDERR_TAIL_LINES = 80


def _pyenv_shim() -> str:
    return (
        "if command -v pyenv >/dev/null 2>&1; then "
        'function pyenv() { case "$1" in rehash|shell) return 0;; esac; command pyenv "$@"; }; '
        "fi; "
    )


_THREAD_LIMIT_VARS = (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
    "OMP_NUM_THREADS",
)


def _apply_thread_limits(env: Dict[str, str], omp_threads: int = 1) -> None:
    """Forcibly set all BLAS/OpenMP thread-count vars in env.

    Called after every env construction so that login-script resets
    (conda init, thisbdm, etc.) cannot survive into subprocesses.
    OMP_NUM_THREADS is set separately per ABM replicate via omp_threads.
    """
    for var in _THREAD_LIMIT_VARS:
        env[var] = str(omp_threads)


def bootstrap_bdm(bdm_env_script: Path) -> Dict[str, str]:
    env = os.environ.copy()
    if env.get("BDMSYS"):
        _apply_thread_limits(env)
        return env
    cmd = (
        f"{_pyenv_shim()}"
        f"export BDM_THISBDM_SILENT=true; "
        f"source {shlex.quote(str(bdm_env_script.resolve()))} && "
        'python3 -c "import os,json;print(json.dumps(dict(os.environ)))"'
    )
    proc = subprocess.run(["bash", "--norc", "--noprofile", "-c", cmd], capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Failed to load BioDynaMo from {bdm_env_script}:\n"
            f"{proc.stderr.strip() or proc.stdout.strip()}"
        )
    env.update(json.loads(proc.stdout))
    _apply_thread_limits(env)  # overwrite any login-script resets
    return env


def _drain_stderr_tail(proc: subprocess.Popen[str], tail: Deque[str]) -> None:
    if proc.stderr is None:
        return
    try:
        for line in proc.stderr:
            tail.append(line)
    finally:
        proc.stderr.close()


def _terminate_process_tree(proc: Optional[subprocess.Popen]) -> None:
    if proc is None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    except Exception:
        try:
            proc.terminate()
        except Exception:
            pass
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass
        try:
            proc.wait(timeout=2)
        except Exception:
            pass


def _write_run_summary(rep_dir: Path, payload: Dict[str, object]) -> None:
    (rep_dir / "run_summary.json").write_text(
        json.dumps(payload, indent=2, default=str),
        encoding="utf-8",
    )


def _write_run_error_log(
    rep_dir: Path,
    *,
    returncode: int,
    cmd: Sequence[str],
    stderr_tail: Sequence[str],
    note: str,
) -> None:
    body = [
        f"returncode={returncode}",
        f"command={' '.join(cmd)}",
        f"note={note}",
        "--- stderr (last {0} lines) ---".format(len(stderr_tail)),
        *stderr_tail,
        "---",
        "Re-run a single eval with: python3 scripts/calibrate_control_growth.py ... --verbose-abm",
    ]
    (rep_dir / "run_error.log").write_text("\n".join(body), encoding="utf-8")


def _remove_heavy_logs(rep_dir: Path) -> None:
    for name in ("run.log",):
        path = rep_dir / name
        if path.exists():
            path.unlink()


def render_progress(label: str, fraction: float, *, sim_h: float | None = None, end_h: float | None = None) -> None:
    fraction = max(0.0, min(1.0, fraction))
    filled = int(round(_PROGRESS_WIDTH * fraction))
    bar = "=" * filled + (">" if filled < _PROGRESS_WIDTH else "") + " " * (max(_PROGRESS_WIDTH - filled - 1, 0))
    note = f"  {sim_h:.2f}/{end_h:.1f} h" if sim_h is not None and end_h and end_h > 0 else ""
    print(f"\r  {label} [{bar[:_PROGRESS_WIDTH]}] {fraction * 100:5.1f}%{note}   ", end="", flush=True)


def _short_param(name: str) -> str:
    return name.replace("cancer_cell/", "")


def _format_param_delta(
    fit: FitParameterization,
    x: np.ndarray,
    reference_x: np.ndarray,
    reference_label: str,
) -> List[str]:
    lines: List[str] = []
    for index, name in enumerate(fit.names):
        value = float(x[index])
        ref = float(reference_x[index])
        delta = value - ref
        if abs(delta) < 1.0e-12:
            change = "(unchanged)"
        else:
            change = f"({delta:+.4g} vs {reference_label})"
        lines.append(f"    {_short_param(name):32s} {ref:8.4g} -> {value:8.4g}  {change}")
    return lines


def _compact_param_delta(fit: FitParameterization, x: np.ndarray, reference_x: np.ndarray) -> str:
    parts: List[str] = []
    for index, name in enumerate(fit.names):
        value = float(x[index])
        ref = float(reference_x[index])
        if abs(value - ref) < 1.0e-12:
            continue
        parts.append(f"{_short_param(name)} {ref:.3g}->{value:.3g}")
    return ", ".join(parts) if parts else "unchanged vs reference"


def _format_viability_vs_excel(
    target_times: np.ndarray,
    target_curve: pd.DataFrame,
    preds: np.ndarray,
) -> List[str]:
    excel = {
        float(row["time_h"]): float(row["target_viability_pct"])
        for _, row in target_curve.iterrows()
    }
    lines = ["  viability (ABM vs Excel, % of t0):"]
    for time_h in target_times:
        pred = prediction_at_time(target_times, preds, float(time_h))
        target = excel.get(float(time_h), float("nan"))
        if not np.isfinite(pred) or not np.isfinite(target):
            continue
        lines.append(
            f"    {time_h:4.0f} h:  ABM {pred:7.1f}%   Excel {target:7.1f}%   Δ {pred - target:+7.1f}"
        )
    return lines


class ControlGrowthObjective:
    def __init__(
        self,
        *,
        abm_binary: Path,
        template_input: Path,
        run_root: Path,
        target_curve: pd.DataFrame,
        fit: Sequence[str],
        time_step_h: float,
        omp_threads: int,
        seed: int,
        replicates: int,
        replicate_seed_step: int,
        bdm_env_script: Path,
        overgrowth_factor: float,
        replicate_timeout_s: int,
        population_band_std_scale: float,
        truncation_penalty_scale: float,
        early_extinction_penalty_scale: float,
        late_extinction_penalty_scale: float,
        late_viability_floor_pct: float,
        verbose_abm: bool = False,
        save_abm_full_log: bool = False,
        max_evals: int = 0,
    ) -> None:
        self.abm_binary = abm_binary.resolve()
        self.template_input = template_input.resolve()
        self.run_root = run_root.resolve()
        self.target_curve = target_curve.sort_values("time_h").reset_index(drop=True)
        self.target_times = self.target_curve["time_h"].to_numpy(dtype=float)
        self.required_end_h = float(np.max(self.target_times))
        positive = self.target_times[self.target_times > 0.0]
        self.first_alive_h = float(positive[0]) if len(positive) else self.required_end_h
        self.fit = select_fit_parameters(fit)
        self.time_step_h = time_step_h
        self.omp_threads = omp_threads
        self.seed = seed
        self.replicates = replicates
        self.replicate_seed_step = replicate_seed_step
        self.replicate_seeds = tuple(seed + i * replicate_seed_step for i in range(replicates))
        self.bdm_env = bootstrap_bdm(bdm_env_script)
        self.verbose_abm = verbose_abm
        self.save_abm_full_log = save_abm_full_log
        self.base_rows = read_rows(self.template_input)
        self.initial_population = int(
            round(get_float_param(self.base_rows, "cancer_cell/initial_population", 1000.0))
        )
        self.template_x = template_point(self.base_rows, self.fit)
        bound_times, floor_counts, ceiling_counts = population_bounds(
            self.target_curve,
            self.initial_population,
            population_band_std_scale,
        )
        self.floor_times = bound_times
        self.floor_counts = floor_counts
        self.ceiling_times = bound_times
        self.ceiling_counts = ceiling_counts
        max_target_viability = float(
            np.max(self.target_curve["target_viability_pct"].to_numpy(dtype=float))
        ) if "target_viability_pct" in self.target_curve.columns else 100.0
        auto_limit = int(math.ceil(
            self.initial_population * (max_target_viability / 100.0) * overgrowth_factor
        ))
        self.total_cell_limit = max(
            auto_limit,
            self.initial_population * 2,
            int(math.ceil(float(np.max(ceiling_counts)))) if len(ceiling_counts) else auto_limit,
        )
        print(
            f"  [cell limit] initial={self.initial_population}  "
            f"max_target_viability={max_target_viability:.1f}%  "
            f"overgrowth_factor={overgrowth_factor:.1f}x  "
            f"-> total_cell_limit={self.total_cell_limit}",
            flush=True,
        )
        self.replicate_timeout_s = max(int(replicate_timeout_s), 0)
        self.truncation_penalty_scale = truncation_penalty_scale
        self.early_extinction_penalty_scale = early_extinction_penalty_scale
        self.late_extinction_penalty_scale = late_extinction_penalty_scale
        self.late_viability_floor_pct = late_viability_floor_pct
        self.eval_id = 0
        self.max_evals = max_evals
        self.current_method = "template"
        self.method_eval_counts: Dict[str, int] = {"template": 0}
        self.best_score = float("inf")
        self.best_x: Optional[np.ndarray] = None
        self.best_pred: Optional[np.ndarray] = None
        self.best_details: Optional[Dict[str, object]] = None
        self._last_eval_x: Optional[np.ndarray] = None
        self.cache: Dict[Tuple[float, ...], Tuple[float, np.ndarray, Dict[str, object]]] = {}
        self.run_root.mkdir(parents=True, exist_ok=True)

    def set_method(self, method: str) -> None:
        self.current_method = method
        self.method_eval_counts.setdefault(method, 0)

    def _run_replicate(self, x: np.ndarray, eval_dir: Path, eval_index: int, rep_index: int, rep_seed: int) -> Tuple[float, np.ndarray, Dict[str, object]]:
        rep_dir = eval_dir / f"rep_{rep_index:02d}"
        if rep_dir.exists():
            import shutil

            shutil.rmtree(rep_dir)
        rep_dir.mkdir(parents=True)
        rows = [row.copy() for row in self.base_rows]
        bounded = clip_to_bounds(np.asarray(x, dtype=float), self.fit.bounds)
        _INT_FIT_PARAMS = {"cancer_cell/can_divide/time_window"}
        for index, name in enumerate(self.fit.names):
            if name in _INT_FIT_PARAMS:
                set_param(rows, name, str(int(round(float(bounded[index])))), "int")
            else:
                set_param(rows, name, f"{float(bounded[index]):.8g}")
        total_steps = int(round(self.required_end_h / self.time_step_h))
        set_param(rows, "output_directory", str(rep_dir / "results"))
        set_param(rows, "number_of_time_steps", str(total_steps), "int")
        set_param(rows, "time_step", f"{self.time_step_h:g}")
        set_param(rows, "statistics_interval", "1", "int")
        set_param(rows, "visualization_interval", str(total_steps + 1), "int")
        set_param(rows, "cancer_cell/initial_population", str(self.initial_population), "int")
        set_param(rows, "clean_output_directory", "true", "bool")
        set_param(rows, "simulation/early_stop_on_total_cells_exceeded", "true", "bool")
        set_param(rows, "simulation/total_cell_limit", str(self.total_cell_limit), "int")
        input_path = rep_dir / "input.csv"
        write_rows(input_path, rows)

        stats_path = rep_dir / "results" / "stats.csv"
        env = self.bdm_env.copy()
        _apply_thread_limits(env, omp_threads=1)  # keep Python-side BLAS single-threaded
        env["OMP_NUM_THREADS"] = str(self.omp_threads)  # ABM binary gets its own OMP budget
        cmd = [str(self.abm_binary), str(input_path), str(rep_seed)]
        label = f"eval {eval_index:04d} rep{rep_index:02d}"
        started = time.time()
        returncode = -1
        timed_out = False
        overgrowth_killed = False
        stderr_tail: Deque[str] = collections.deque(maxlen=_STDERR_TAIL_LINES)
        stderr_thread: Optional[threading.Thread] = None
        proc: Optional[subprocess.Popen] = None
        try:
            if self.verbose_abm:
                proc = subprocess.Popen(cmd, cwd=self.template_input.parent, env=env, start_new_session=True)
                while proc.poll() is None:
                    if self.replicate_timeout_s > 0 and (time.time() - started) > self.replicate_timeout_s:
                        timed_out = True
                        _terminate_process_tree(proc)
                        break
                    _ch, _nc = probe_stats_cells(stats_path)
                    if _nc > self.total_cell_limit and _ch < self.required_end_h * 0.8:
                        overgrowth_killed = True
                        _terminate_process_tree(proc)
                        print(f"\n  !! early overgrowth kill at {_ch:.1f}h: N_cells={_nc} > {self.total_cell_limit}", flush=True)
                        break
                    time.sleep(0.2)
                returncode = int(proc.returncode or -1)
            elif self.save_abm_full_log:
                log_path = rep_dir / "run.log"
                with log_path.open("w", encoding="utf-8") as log:
                    proc = subprocess.Popen(
                        cmd,
                        cwd=self.template_input.parent,
                        env=env,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        start_new_session=True,
                    )
                    while proc.poll() is None:
                        if self.replicate_timeout_s > 0 and (time.time() - started) > self.replicate_timeout_s:
                            timed_out = True
                            _terminate_process_tree(proc)
                            break
                        _ch, _nc = probe_stats_cells(stats_path)
                        if _nc > self.total_cell_limit and _ch < self.required_end_h * 0.8:
                            overgrowth_killed = True
                            _terminate_process_tree(proc)
                            print(f"\n  !! early overgrowth kill at {_ch:.1f}h: N_cells={_nc} > {self.total_cell_limit}", flush=True)
                            break
                        sim_h, frac = probe_stats_progress(stats_path, self.required_end_h)
                        render_progress(label, frac, sim_h=sim_h, end_h=self.required_end_h)
                        time.sleep(0.2)
                    returncode = int(proc.returncode or -1)
            else:
                proc = subprocess.Popen(
                    cmd,
                    cwd=self.template_input.parent,
                    env=env,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                    text=True,
                    bufsize=1,
                    start_new_session=True,
                )
                stderr_thread = threading.Thread(
                    target=_drain_stderr_tail,
                    args=(proc, stderr_tail),
                    daemon=True,
                )
                stderr_thread.start()
                while proc.poll() is None:
                    if self.replicate_timeout_s > 0 and (time.time() - started) > self.replicate_timeout_s:
                        timed_out = True
                        _terminate_process_tree(proc)
                        break
                    _ch, _nc = probe_stats_cells(stats_path)
                    if _nc > self.total_cell_limit and _ch < self.required_end_h * 0.8:
                        overgrowth_killed = True
                        _terminate_process_tree(proc)
                        print(f"\n  !! early overgrowth kill at {_ch:.1f}h: N_cells={_nc} > {self.total_cell_limit}", flush=True)
                        break
                    sim_h, frac = probe_stats_progress(stats_path, self.required_end_h)
                    render_progress(label, frac, sim_h=sim_h, end_h=self.required_end_h)
                    time.sleep(0.2)
                returncode = int(proc.returncode or -1)
                if stderr_thread is not None:
                    stderr_thread.join(timeout=5.0)
            if not self.verbose_abm:
                sim_h, frac = probe_stats_progress(stats_path, self.required_end_h)
                render_progress(
                    label,
                    1.0 if returncode == 0 and sim_h >= self.required_end_h - 1.0e-6 else frac,
                    sim_h=sim_h,
                    end_h=self.required_end_h,
                )
                print(flush=True)
        except KeyboardInterrupt:
            _terminate_process_tree(proc)
            if stderr_thread is not None:
                stderr_thread.join(timeout=2.0)
            raise
        except Exception:
            _terminate_process_tree(proc)
            _write_run_error_log(
                rep_dir,
                returncode=-1,
                cmd=cmd,
                stderr_tail=list(stderr_tail),
                note="subprocess_exception",
            )
            nan = np.full(len(self.target_times), np.nan)
            return 1.0e6, nan, {"termination_state": "subprocess_exception", "returncode": -1}

        if not stats_path.exists():
            _terminate_process_tree(proc)
            _write_run_summary(
                rep_dir,
                {
                    "returncode": returncode,
                    "seed": rep_seed,
                    "elapsed_s": time.time() - started,
                    "termination_state": "timeout" if timed_out else "missing_stats",
                },
            )
            _write_run_error_log(
                rep_dir,
                returncode=returncode,
                cmd=cmd,
                stderr_tail=list(stderr_tail),
                note="timeout" if timed_out else "missing_stats",
            )
            compact_output_dir(rep_dir / "results")
            nan = np.full(len(self.target_times), np.nan)
            return 1.0e6, nan, {
                "termination_state": "timeout" if timed_out else "missing_stats",
                "returncode": returncode,
            }

        stats = pd.read_csv(stats_path, skipinitialspace=True)
        times = stats["current_time"].to_numpy(dtype=float)
        viable = viable_counts(stats)
        total = stats["N_cells"].to_numpy(dtype=float) if "N_cells" in stats.columns else viable
        score, preds, details = evaluate_growth_trace(
            times,
            viable,
            total,
            self.target_curve,
            self.target_times,
            self.required_end_h,
            self.first_alive_h,
            self.initial_population,
            self.floor_times,
            self.floor_counts,
            self.ceiling_times,
            self.ceiling_counts,
            truncation_scale=self.truncation_penalty_scale,
            early_extinction_scale=self.early_extinction_penalty_scale,
            late_extinction_scale=self.late_extinction_penalty_scale,
            late_viability_floor_pct=self.late_viability_floor_pct,
        )
        details["returncode"] = returncode
        details["seed"] = rep_seed
        details["elapsed_s"] = time.time() - started
        details["timed_out"] = timed_out
        details["overgrowth_killed"] = overgrowth_killed
        summary = {
            "returncode": returncode,
            "seed": rep_seed,
            "elapsed_s": details["elapsed_s"],
            "termination_state": details.get("termination_state"),
            "score": score,
            "final_time_h": details.get("final_time_h"),
            "predicted_viability_pct": preds.tolist(),
        }
        _write_run_summary(rep_dir, summary)
        failed = returncode != 0 or bool(details.get("invalid_metrics"))
        if failed:
            _write_run_error_log(
                rep_dir,
                returncode=returncode,
                cmd=cmd,
                stderr_tail=list(stderr_tail),
                note=str(details.get("termination_state", "failed")),
            )
        elif not self.save_abm_full_log:
            _remove_heavy_logs(rep_dir)
        compact_output_dir(rep_dir / "results")
        return score, preds, details

    def _print_eval_banner(
        self,
        *,
        method: str,
        method_eval_id: int,
        eval_index: int,
        x: np.ndarray,
        cached: bool = False,
    ) -> None:
        tag = f"[{method} #{method_eval_id:03d}] eval {eval_index:04d}"
        if cached:
            print(f"\n{tag}  (cached — skipping ABM reruns)", flush=True)
            return
        print(f"\n{tag}", flush=True)
        print("  parameters vs template:", flush=True)
        for line in _format_param_delta(self.fit, x, self.template_x, "template"):
            print(line, flush=True)
        if self._last_eval_x is not None:
            print("  parameters vs previous eval:", flush=True)
            for line in _format_param_delta(self.fit, x, self._last_eval_x, "previous"):
                print(line, flush=True)
        compact = _compact_param_delta(self.fit, x, self.template_x)
        print(f"  Δ summary: {compact}", flush=True)

    def _print_eval_result(
        self,
        *,
        method: str,
        method_eval_id: int,
        eval_index: int,
        x: np.ndarray,
        mean_score: float,
        score_std: float,
        mean_pred: np.ndarray,
        elapsed_s: float,
        cached: bool = False,
        is_new_best: bool = False,
    ) -> None:
        best_tag = "  ★ new best" if is_new_best and not cached else ""
        cache_tag = "  [cached]" if cached else ""
        print(
            f"  → score {mean_score:.4f} ± {score_std:.4f}  |  "
            f"best {min(self.best_score, mean_score):.4f}  |  {elapsed_s:.0f}s{best_tag}{cache_tag}",
            flush=True,
        )
        for line in _format_viability_vs_excel(self.target_times, self.target_curve, mean_pred):
            print(line, flush=True)

    def simulate(self, x: np.ndarray) -> Tuple[float, np.ndarray, Dict[str, object]]:
        bounded_x = clip_to_bounds(np.asarray(x, dtype=float), self.fit.bounds)
        if self.max_evals > 0 and self.eval_id >= self.max_evals:
            dummy_pred = np.full(len(self.target_times), float("nan"))
            return 1e9, dummy_pred, {"budget_exhausted": True, "eval_id": self.eval_id}
        key = tuple(float(v) for v in bounded_x)
        eval_index = self.eval_id
        method = self.current_method
        method_eval_id = self.method_eval_counts.get(method, 0)

        if key in self.cache:
            mean_score, mean_pred, details = self.cache[key]
            if mean_score < self.best_score:
                self.best_score = mean_score
                self.best_x = bounded_x.copy()
                self.best_pred = mean_pred.copy()
                self.best_details = dict(details)
            self._print_eval_banner(
                method=method,
                method_eval_id=method_eval_id,
                eval_index=eval_index,
                x=bounded_x,
                cached=True,
            )
            self._print_eval_result(
                method=method,
                method_eval_id=method_eval_id,
                eval_index=eval_index,
                x=bounded_x,
                mean_score=mean_score,
                score_std=float(details.get("score_std", 0.0)),
                mean_pred=mean_pred,
                elapsed_s=0.0,
                cached=True,
            )
            return self.cache[key]

        self.method_eval_counts[method] = method_eval_id + 1
        self._print_eval_banner(
            method=method,
            method_eval_id=method_eval_id,
            eval_index=eval_index,
            x=bounded_x,
        )
        eval_dir = self.run_root / f"eval_{eval_index:04d}"
        if eval_dir.exists():
            import shutil

            shutil.rmtree(eval_dir)
        eval_dir.mkdir(parents=True)

        started = time.time()
        rep_scores: List[float] = []
        rep_preds: List[np.ndarray] = []
        rep_details: List[Dict[str, object]] = []
        for rep_index, rep_seed in enumerate(self.replicate_seeds):
            score, preds, details = self._run_replicate(bounded_x, eval_dir, eval_index, rep_index, rep_seed)
            rep_scores.append(score)
            rep_preds.append(preds)
            rep_details.append(details)

        score_arr = np.asarray(rep_scores, dtype=float)
        pred_arr = np.stack(rep_preds, axis=0)
        mean_score = float(np.mean(score_arr))
        mean_pred = np.nanmean(pred_arr, axis=0)
        score_std = float(np.std(score_arr)) if len(score_arr) > 1 else 0.0
        is_new_best = mean_score < self.best_score - 1.0e-12
        if is_new_best:
            self.best_score = mean_score
            self.best_x = bounded_x.copy()
            self.best_pred = mean_pred.copy()
            self.best_details = None
        else:
            self.best_score = min(self.best_score, mean_score)
        param_snapshot = {self.fit.names[i]: float(bounded_x[i]) for i in range(len(bounded_x))}
        param_vs_template = {
            self.fit.names[i]: {
                "template": float(self.template_x[i]),
                "value": float(bounded_x[i]),
                "delta": float(bounded_x[i] - self.template_x[i]),
            }
            for i in range(len(bounded_x))
        }

        primary = rep_details[0] if rep_details else {}
        details = {
            **primary,
            "method": method,
            "method_eval_id": method_eval_id,
            "eval_id": eval_index,
            "parameters": param_snapshot,
            "parameters_vs_template": param_vs_template,
            "replicate_count": self.replicates,
            "replicate_scores": score_arr.tolist(),
            "score_std": score_std,
            "best_score_global": self.best_score,
            "predicted_viability_pct": mean_pred.tolist(),
            "replicate_details": rep_details,
        }
        self.cache[key] = (mean_score, mean_pred, details)
        self.eval_id += 1
        self._last_eval_x = bounded_x.copy()

        self._print_eval_result(
            method=method,
            method_eval_id=method_eval_id,
            eval_index=eval_index,
            x=bounded_x,
            mean_score=mean_score,
            score_std=score_std,
            mean_pred=mean_pred,
            elapsed_s=time.time() - started,
            is_new_best=is_new_best,
        )
        return self.cache[key]

    def __call__(self, x: np.ndarray) -> float:
        return self.simulate(x)[0]


def run_bayesian_optimization(
    objective: ControlGrowthObjective,
    bounds: Sequence[Tuple[float, float]],
    seed: int,
    n_calls: int = 100,
    n_initial_points: int = 20,
) -> Tuple[np.ndarray, float, dict]:
    """Bayesian Optimization with Gaussian Process surrogate for expensive ABM.
    
    Handles noisy objectives via the 'gaussian' noise option.
    Returns: (best_x, best_score, result_dict)
    """
    if not SKOPT_AVAILABLE:
        raise RuntimeError("scikit-optimize not installed. Run: pip install scikit-optimize")
    
    # Use module-level imported functions
    gp_minimize = _gp_minimize
    Real = _Real
    use_named_args = _use_named_args
    
    # Define search space
    space = [Real(lo, hi, name=name) for name, (lo, hi) in zip(objective.fit.names, bounds)]
    
    # Wrap objective to return scalar (skopt expects this)
    @use_named_args(space)
    def skopt_objective(**kwargs):
        x = np.array([kwargs[name] for name in objective.fit.names])
        score, _, _ = objective.simulate(x)
        return score
    
    print(f"\n=== Bayesian Optimization (GP surrogate) ===", flush=True)
    print(f"  Total budget: {n_calls} evaluations", flush=True)
    print(f"  Random initialization: {n_initial_points} points", flush=True)
    print(f"  Acquisition: Expected Improvement (EI)", flush=True)
    print(f"  Noise handling: Gaussian", flush=True)

    gp_kwargs = {
        "n_calls": n_calls,
        "n_initial_points": n_initial_points,
        "noise": "gaussian",  # Critical: handles stochastic ABM replicates
        "acq_func": "EI",     # Expected improvement
        "random_state": seed,
        "verbose": True,
        "n_jobs": 1,           # ABM evaluations are sequential in this script
    }
    if threadpool_limits is None:
        result = gp_minimize(skopt_objective, space, **gp_kwargs)
    else:
        with threadpool_limits(limits=1):
            result = gp_minimize(skopt_objective, space, **gp_kwargs)
    
    best_x = clip_to_bounds(np.asarray(result.x, dtype=float), bounds)
    best_score = float(result.fun)
    
    result_dict = {
        "convergence_iters": len(result.func_vals),
        "models": len(result.models) if hasattr(result, 'models') else 0,
        "space": [str(s) for s in space],
    }
    
    return best_x, best_score, result_dict


def run_optimizers(
    objective: ControlGrowthObjective,
    methods: Sequence[str],
    seed: int,
    *,
    da_maxiter: int,
    nm_maxiter: int,
    max_evals: int,
    bo_calls: int = 100,
    bo_init: int = 20,
) -> Tuple[Dict[str, dict], np.ndarray]:
    results: Dict[str, dict] = {}
    bounds = list(objective.fit.bounds)
    template_x = clip_to_bounds(objective.template_x.copy(), bounds)
    print("\n=== Template baseline (input_control.csv values) ===", flush=True)
    for index, name in enumerate(objective.fit.names):
        print(f"  {_short_param(name):32s} {template_x[index]:8.4g}", flush=True)

    objective.set_method("template")
    template_score, template_pred, template_details = objective.simulate(template_x)
    results["template"] = {
        "score": template_score,
        "parameters": {objective.fit.names[i]: float(template_x[i]) for i in range(len(template_x))},
        "predicted_viability_pct": template_pred.tolist(),
        "details": template_details,
    }
    incumbent_x = template_x.copy()
    incumbent_score = template_score

    for method in methods:
        objective.set_method(method)
        print(f"\n=== Optimizer: {method} ===", flush=True)
        if method == "bayesian":
            if not SKOPT_AVAILABLE:
                print("  WARNING: scikit-optimize not available, skipping bayesian", flush=True)
                continue
            x, score, details = run_bayesian_optimization(
                objective, bounds, seed, n_calls=bo_calls, n_initial_points=bo_init
            )
            bo_key = tuple(float(v) for v in x)
            if bo_key in objective.cache:
                _bo_score, bo_pred, _bo_details = objective.cache[bo_key]
            else:
                bo_pred = objective.best_pred
            results[method] = {
                "score": score,
                "parameters": {objective.fit.names[i]: float(x[i]) for i in range(len(x))},
                "predicted_viability_pct": bo_pred.tolist() if bo_pred is not None else [],
                "details": details,
            }
            if score <= incumbent_score:
                incumbent_x = x.copy()
                incumbent_score = score
            continue
        elif method == "dual_annealing":
            res = dual_annealing(objective, bounds=bounds, seed=seed + 11, maxiter=da_maxiter, maxfun=max_evals * 10)
            x = clip_to_bounds(np.asarray(res.x, dtype=float), bounds)
        elif method == "powell":
            res = minimize(
                objective,
                x0=incumbent_x,
                method="Powell",
                bounds=bounds,
                options={"maxiter": nm_maxiter, "maxfev": max_evals * 10},
            )
            x = clip_to_bounds(np.asarray(res.x, dtype=float), bounds)
        elif method == "nelder_mead":
            res = minimize(
                objective,
                x0=incumbent_x,
                method="Nelder-Mead",
                bounds=bounds,
                options={"maxiter": nm_maxiter, "maxfev": max_evals * 10},
            )
            x = clip_to_bounds(np.asarray(res.x, dtype=float), bounds)
        else:
            raise ValueError(f"Unsupported method: {method}")

        key = tuple(float(v) for v in x)
        if key in objective.cache:
            score, pred, details = objective.cache[key]
        elif objective.best_x is not None:
            best_key = tuple(float(v) for v in objective.best_x)
            if best_key in objective.cache:
                score, pred, details = objective.cache[best_key]
                x = objective.best_x.copy()
            else:
                score = objective.best_score
                pred = objective.best_pred if objective.best_pred is not None else np.full(len(objective.target_times), float("nan"))
                details = objective.best_details or {}
                x = objective.best_x.copy()
        else:
            score = objective.best_score
            pred = objective.best_pred if objective.best_pred is not None else np.full(len(objective.target_times), float("nan"))
            details = objective.best_details or {}
        results[method] = {
            "score": score,
            "parameters": {objective.fit.names[i]: float(x[i]) for i in range(len(x))},
            "predicted_viability_pct": pred.tolist(),
            "details": details,
        }
        if score <= incumbent_score:
            incumbent_x = x.copy()
            incumbent_score = score
    return results, incumbent_x


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--excel", type=Path, required=True)
    parser.add_argument("--cell-line", default="EGI1")
    parser.add_argument("--template-input", type=Path, default=Path("input_control.csv"))
    parser.add_argument("--abm-binary", type=Path, default=Path("../../build/ABM4bio"))
    parser.add_argument("--bdm-env-script", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("calibration_outputs/EGI1_control_0to72h"))
    parser.add_argument("--fit-params", nargs="+", default=list(CONTROL_GROWTH_PARAMS.names[:4]))
    parser.add_argument("--methods", nargs="+", default=["bayesian", "powell"])
    parser.add_argument("--bo-calls", type=int, default=80, help="Bayesian opt: total function evaluations")
    parser.add_argument("--bo-init", type=int, default=20, help="Bayesian opt: random initialization points")
    parser.add_argument("--time-step-h", type=float, default=0.01)
    parser.add_argument("--omp-threads", type=int, default=2)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--replicate-seed-step", type=int, default=1)
    parser.add_argument("--da-maxiter", type=int, default=3)
    parser.add_argument("--nm-maxiter", type=int, default=20)
    parser.add_argument("--max-evaluations", type=int, default=200)
    parser.add_argument(
        "--overgrowth-factor",
        type=float,
        default=4.0,
        help="Kill replicate when N_cells > initial_population × max_target_viability × this factor (default: 4.0)",
    )
    parser.add_argument("--replicate-timeout-s", type=int, default=900)
    parser.add_argument("--population-band-std-scale", type=float, default=1.0)
    parser.add_argument("--truncation-penalty-scale", type=float, default=12.0)
    parser.add_argument("--early-extinction-penalty-scale", type=float, default=8.0)
    parser.add_argument("--late-extinction-penalty-scale", type=float, default=8.0)
    parser.add_argument("--late-viability-floor-pct", type=float, default=25.0)
    parser.add_argument("--verbose-abm", action="store_true")
    parser.add_argument(
        "--save-abm-full-log",
        action="store_true",
        help="Write full ABM stdout to run.log per replicate (large; default is summary only).",
    )
    args = parser.parse_args()

    if args.seed <= 0:
        raise SystemExit("--seed must be > 0")
    if args.replicates < 1:
        raise SystemExit("--replicates must be >= 1")

    df = parse_workbook(args.excel)
    target = build_control_target(df, args.cell_line, normalization="t0")
    print(f"Target: {args.cell_line} untreated control ({len(target)} time points, 0–72 h, t0 norm)")
    for _, row in target.iterrows():
        print(f"  t={row['time_h']:g}h  target={row['target_viability_pct']:.1f}% ± {row['target_sd_pct']:.1f}")

    objective = ControlGrowthObjective(
        abm_binary=args.abm_binary,
        template_input=args.template_input,
        run_root=args.out_dir / "abm_runs",
        target_curve=target,
        fit=args.fit_params,
        time_step_h=args.time_step_h,
        omp_threads=args.omp_threads,
        seed=args.seed,
        replicates=args.replicates,
        replicate_seed_step=args.replicate_seed_step,
        bdm_env_script=args.bdm_env_script,
        overgrowth_factor=args.overgrowth_factor,
        replicate_timeout_s=args.replicate_timeout_s,
        population_band_std_scale=args.population_band_std_scale,
        truncation_penalty_scale=args.truncation_penalty_scale,
        early_extinction_penalty_scale=args.early_extinction_penalty_scale,
        late_extinction_penalty_scale=args.late_extinction_penalty_scale,
        late_viability_floor_pct=args.late_viability_floor_pct,
        verbose_abm=args.verbose_abm,
        save_abm_full_log=args.save_abm_full_log,
        max_evals=args.max_evaluations,
    )

    try:
        results, best_x = run_optimizers(
            objective,
            args.methods,
            args.seed,
            da_maxiter=args.da_maxiter,
            nm_maxiter=args.nm_maxiter,
            max_evals=args.max_evaluations,
            bo_calls=args.bo_calls,
            bo_init=args.bo_init,
        )
    except KeyboardInterrupt:
        print("\nInterrupted by user (Ctrl+C). Stopped cleanly.", flush=True)
        raise SystemExit(130)

    best_score = objective.best_score
    sens_rows = sensitivity_report(objective, best_x, best_score)
    print_sensitivity_report(sens_rows)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "cell_line": args.cell_line,
        "normalization": "t0",
        "target_curve": target.to_dict(orient="records"),
        "results": results,
        "best_parameters": {objective.fit.names[i]: float(best_x[i]) for i in range(len(best_x))},
        "sensitivity": sens_rows,
    }
    out_json = args.out_dir / "control_calibration_results.json"
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nWrote {out_json}")


if __name__ == "__main__":
    main()
