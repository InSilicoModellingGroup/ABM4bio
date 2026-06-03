"""Load EGI1 (and other) viability trajectories from the Gorjet Excel workbook."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}

EXPOSURE_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([\"']|min)?\s*$", re.IGNORECASE)
TIME_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*h\s*$", re.IGNORECASE)


def col_ref_to_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    value = 0
    for ch in letters:
        value = value * 26 + (ord(ch.upper()) - ord("A") + 1)
    return value - 1


def parse_exposure_minutes(label: str) -> float | None:
    if not label:
        return None
    match = EXPOSURE_RE.match(label)
    if not match:
        return None
    value = float(match.group(1))
    unit = (match.group(2) or "").lower()
    # All exposure labels in this workbook are in minutes:
    #   "0", "2", "5" (unitless) → minutes
    #   "2'", "5min"             → minutes
    # There are no second or hour exposure labels in this dataset.
    return value  # always minutes


def parse_time_hours(label: str) -> float | None:
    if not label:
        return None
    match = TIME_RE.match(label)
    return float(match.group(1)) if match else None


def load_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    return [
        "".join((node.text or "") for node in item.findall(".//a:t", NS))
        for item in root.findall("a:si", NS)
    ]


def parse_sheet_cells(
    zf: zipfile.ZipFile, sheet_path: str, shared_strings: list[str]
) -> dict[tuple[int, int], str]:
    root = ET.fromstring(zf.read(sheet_path))
    cells: dict[tuple[int, int], str] = {}
    for row in root.findall(".//a:sheetData/a:row", NS):
        row_index = int(row.attrib["r"]) - 1
        for cell in row.findall("a:c", NS):
            cref = cell.attrib.get("r", "")
            col_index = col_ref_to_index(cref) if cref else 0
            cell_type = cell.attrib.get("t", "")
            value_node = cell.find("a:v", NS)
            value = "" if value_node is None or value_node.text is None else value_node.text
            if cell_type == "s" and value:
                value = shared_strings[int(value)]
            cells[(row_index, col_index)] = value.strip()
    return cells


def parse_workbook(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as zf:
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        relationships = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in relationships}
        shared_strings = load_shared_strings(zf)

        records: list[dict] = []
        for sheet in workbook.findall("a:sheets/a:sheet", NS):
            sheet_name = sheet.attrib["name"]
            if "final" not in sheet_name.lower():
                continue
            rel_id = sheet.attrib[f"{{{NS['r']}}}id"]
            sheet_path = "xl/" + rel_map[rel_id]
            cells = parse_sheet_cells(zf, sheet_path, shared_strings)
            for anchor_col in (0, 8):
                for row in range(200):
                    exposure_label = cells.get((row, anchor_col), "")
                    exposure_min = parse_exposure_minutes(exposure_label)
                    if exposure_min is None:
                        continue
                    for offset in range(1, 5):
                        time_label = cells.get((row + offset, anchor_col), "")
                        time_h = parse_time_hours(time_label)
                        if time_h is None:
                            continue
                        mean_raw = cells.get((row + offset, anchor_col + 5), "")
                        sd_raw = cells.get((row + offset, anchor_col + 6), "")
                        if not mean_raw:
                            continue
                        records.append(
                            {
                                "cell_line": sheet_name.replace(" (FINAL)", "").strip(),
                                "exposure_min": float(exposure_min),
                                "time_h": float(time_h),
                                "mean_signal": float(mean_raw),
                                "sd_signal": float(sd_raw) if sd_raw else np.nan,
                            }
                        )

    df = pd.DataFrame.from_records(records).drop_duplicates(
        subset=["cell_line", "exposure_min", "time_h"], keep="last"
    )
    df.sort_values(["cell_line", "exposure_min", "time_h"], inplace=True)
    return df.reset_index(drop=True)


def build_control_target(
    df_all: pd.DataFrame,
    cell_line: str,
    *,
    normalization: str = "t0",
) -> pd.DataFrame:
    """Untreated control curve (0 min CAP) for 0–72 h calibration."""
    df = df_all[df_all["cell_line"].str.lower() == cell_line.lower()].copy()
    if df.empty:
        valid = sorted(df_all["cell_line"].unique())
        raise ValueError(f'Cell line "{cell_line}" not found. Options: {valid}')

    available = np.sort(df["exposure_min"].unique())
    control_durations = available[np.isclose(available, 0.0, atol=1.0e-12)]
    if len(control_durations) == 0:
        raise ValueError(f"No untreated control (0 min exposure) data for {cell_line}.")
    used_duration = float(control_durations[0])
    curve = df[np.isclose(df["exposure_min"], used_duration)].copy().sort_values("time_h")
    if curve.empty:
        raise ValueError(f"No control time points for {cell_line}.")

    if normalization == "t0":
        baseline = float(curve.loc[np.isclose(curve["time_h"], 0.0), "mean_signal"].iloc[0])
        curve["target_viability_pct"] = 100.0 * curve["mean_signal"] / baseline
        curve["target_sd_pct"] = 100.0 * curve["sd_signal"] / baseline
    else:
        controls = (
            df_all[df_all["exposure_min"] == 0.0][["cell_line", "time_h", "mean_signal", "sd_signal"]]
            .rename(columns={"mean_signal": "ctrl_mean", "sd_signal": "ctrl_sd"})
        )
        curve = curve.merge(controls, on=["cell_line", "time_h"], how="left")
        curve["target_viability_pct"] = 100.0 * curve["mean_signal"] / curve["ctrl_mean"]
        curve["target_sd_pct"] = 100.0 * np.sqrt(
            (curve["sd_signal"] / curve["ctrl_mean"]) ** 2
            + ((curve["mean_signal"] * curve["ctrl_sd"]) / (curve["ctrl_mean"] ** 2)) ** 2
        )

    curve["target_sd_pct"] = (
        curve["target_sd_pct"].replace([np.inf, -np.inf], np.nan).fillna(5.0).clip(lower=2.0)
    )
    curve["time_days"] = curve["time_h"] / 24.0
    return curve.reset_index(drop=True)
