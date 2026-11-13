#!/usr/bin/env python3
"""Compare regenerated outputs against selected numerical claims in final_2.pdf."""
import argparse
import csv
import json
import platform
from pathlib import Path


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanisms", type=Path)
    parser.add_argument("--h1", type=Path)
    parser.add_argument("--phase29", type=Path)
    parser.add_argument("--phase28", type=Path)
    parser.add_argument("--phase28-history", type=Path)
    parser.add_argument("--forecast", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("--out must name a fresh directory")
    args.out.mkdir(parents=True)
    rows = []

    def add(claim, expected, actual, tolerance, level, source):
        if actual is None:
            status, diff, relative = "BLOCKED", None, None
        else:
            diff = abs(actual - expected)
            relative = diff / abs(expected) if expected else None
            status = "PASS" if diff <= tolerance else "FAIL"
        rows.append(dict(claim=claim, expected=expected, actual=actual,
                         absolute_difference=diff, relative_difference=relative,
                         tolerance=tolerance, status=status,
                         reproduction_level=level, source=source))

    def locate(root, filename):
        path = root / filename if root else None
        return path if path and path.is_file() else None

    m = args.mechanisms
    integrity = locate(m, "model_integrity.csv")
    if integrity:
        for item in read_csv(integrity):
            add("M2_REPLAY_" + item["component"], 0,
                float(item["max_absolute_error"]), 1e-6, 2, str(integrity))
    else:
        add("M2_REPLAY", 0, None, 1e-6, 2, "model_integrity.csv")
    table33 = {
        2013: (-25.431, -31.243, -1.165, -0.417),
        2014: (-23.614, -29.685, -1.301, -0.477),
        2015: (-28.241, -34.128, -2.483, -0.879),
        2016: (-24.772, -29.121, -5.389, -2.144),
        2017: (-6.673, -8.497, -2.559, -0.964),
        2018: (-7.584, -9.511, -2.472, -0.925),
        2019: (-9.784, -11.986, -2.495, -0.928),
        2021: (11.102, 7.020, -4.683, None),
    }
    for policy in ("raw", "blend"):
        h2_path = locate(m, f"h2_seasons_{policy}.csv")
        h3_path = locate(m, f"h3_seasons_{policy}.csv")
        finite_path = locate(m, f"h3_finite_{policy}.csv")
        h2 = {int(r["year"]): r for r in read_csv(h2_path)} if h2_path else {}
        h3 = {int(r["year"]): r for r in read_csv(h3_path)} if h3_path else {}
        finite = {int(r["year"]): r for r in read_csv(finite_path)} if finite_path else {}
        add(f"H2_ELIGIBLE_{policy}", 8,
            sum(r["eligible"] == "TRUE" for r in h2.values()) if h2 else None,
            0, 2, str(h2_path or "missing"))
        add(f"H3_LOCAL_ELIGIBLE_{policy}", 8,
            sum(r["eligible"] == "TRUE" for r in h3.values()) if h3 else None,
            0, 2, str(h3_path or "missing"))
        add(f"H3_FINITE_ELIGIBLE_{policy}", 7,
            sum(r["eligible"] == "TRUE" for r in finite.values()) if finite else None,
            0, 2, str(finite_path or "missing"))
        if policy == "raw":
            for year, expected in table33.items():
                for column, value in zip(("E_S", "E_beta"), expected[:2]):
                    add(f"T3.3_{year}_{column}", value,
                        float(h2[year][column]) if year in h2 else None,
                        0.0005, 2, str(h2_path or "missing"))
                add(f"T3.3_{year}_E_I", expected[2],
                    float(h3[year]["E_I"]) if year in h3 else None,
                    0.0005, 2, str(h3_path or "missing"))
                if expected[3] is not None:
                    add(f"T3.3_{year}_FINITE_I", expected[3],
                        float(finite[year]["contrast"]) if year in finite else None,
                        0.0005, 2, str(finite_path or "missing"))
            add("H2_BETA_LARGER", 7,
                sum(abs(float(r["E_beta"])) > abs(float(r["E_S"]))
                    for r in h2.values() if r["eligible"] == "TRUE") if h2 else None,
                0, 2, str(h2_path or "missing"))
    h1_path = locate(args.h1, "h1_parents.csv")
    if h1_path:
        parents = read_csv(h1_path)
        for policy, expected_events, expected_displacements, expected_fraction in (
            ("raw", 80, 72, 0.7083333333333333),
            ("blend", 77, 69, 0.705)):
            group = [r for r in parents if r["policy"] == policy]
            truth = lambda r, key: r[key] == "TRUE"
            events = sum(truth(r, "timing_event") for r in group)
            pure = sum(truth(r, "timing_event") and truth(r, "displacement")
                       and not truth(r, "mixed") and not truth(r, "unresolved")
                       for r in group)
            fractions = []
            for year in range(2013, 2023):
                season = [r for r in group if int(r["season"]) == year]
                n = sum(truth(r, "timing_event") for r in season)
                d = sum(truth(r, "timing_event") and truth(r, "displacement")
                        and not truth(r, "mixed") and not truth(r, "unresolved")
                        for r in season)
                fractions.append(d/n if n else 0)
            add(f"H1_EVENTS_{policy}", expected_events, events, 0, 2, str(h1_path))
            add(f"H1_PURE_{policy}", expected_displacements, pure, 0, 2, str(h1_path))
            add(f"H1_EQUAL_SEASON_{policy}", expected_fraction,
                sum(fractions)/10, 0.00005, 2, str(h1_path))
    else:
        add("H1_EVENTS_raw", 80, None, 0, 2, "h1_parents.csv")
    h1_edges = locate(args.h1, "h1_edges.csv")
    edges = read_csv(h1_edges) if h1_edges else []
    for year, left, right, early, late, gap_left, gap_right in (
        (2017, -0.160600, -0.160244, 24, 11, -4.49, 2.50),
        (2018, -0.169146, -0.168790, 23, 12, -4.30, 4.73),
        (2022, 0.066946, 0.068371, 51, 4, -38.57, 3.02)):
        q = next((r for r in edges if r["policy"] == "raw" and
                  int(r["season"]) == year and
                  r["classification"] == "DOMINANCE_EXCHANGE" and
                  int(r["left_week"]) == early and
                  int(r["right_week"]) == late), None)
        for field, value, tolerance in (
            ("v_left", left, 0.0000005),
            ("v_right", right, 0.0000005),
            ("early_minus_late_left", gap_left, 0.005),
            ("early_minus_late_right", gap_right, 0.005)):
            add(f"T3.2_{year}_{field}", value,
                float(q[field]) if q else None, tolerance, 2,
                str(h1_edges or "h1_edges.csv"))
    phase29 = locate(args.phase29, "fitting_metrics.csv")
    timing29 = locate(args.phase29, "timing_shift.csv")
    metrics29 = read_csv(phase29) if phase29 else []
    for model, mae, rmse in (("SPL_M0", 408.03, 754.24),
                             ("SPL_M1", 319.24, 618.45),
                             ("SPL_M2", 264.74, 533.79),
                             ("A_M2", 360.15, 691.04),
                             ("AS_M2", 314.67, 608.02)):
        q = next((r for r in metrics29 if r["model"] == model and
                  r["period"] == "2013-2023"), None)
        add(f"T3.1_{model}_ROWS", 572, int(q["n"]) if q else None,
            0, 1, str(phase29 or "fitting_metrics.csv"))
        for key, value in (("mae", mae), ("rmse", rmse)):
            add(f"T3.1_{model}_{key.upper()}", value,
                float(q[key]) if q else None, 0.01, 1,
                str(phase29 or "fitting_metrics.csv"))
    shifts = read_csv(timing29) if timing29 else []
    for model, advance in (("SPL_M0", 5), ("SPL_M1", 9),
                           ("SPL_M2", 11), ("A_M2", 9), ("AS_M2", 8)):
        q = next((r for r in shifts if r["model"] == model), None)
        add(f"B.3_{model}_ADVANCE", advance,
            int(q["fitted_advance"]) if q else None, 0, 1,
            str(timing29 or "timing_shift.csv"))
    phase28 = locate(args.phase28, "fit_metrics.csv")
    movement28 = locate(args.phase28, "movement_metrics.csv")
    runtime28 = locate(args.phase28_history or args.phase28, "runtime_metrics.csv")
    inventory28 = locate(args.phase28_history or args.phase28,
                         "parameter_inventory.csv")
    metrics28 = read_csv(phase28) if phase28 else []
    for model, summary, key, expected in (
        ("M2", "M2", "mae", 264.74),
        ("M2", "M2", "rmse", 533.79),
        ("PARK", "PARK_median", "mae", 419.72),
        ("PARK", "PARK_mean", "rmse", 744.74)):
        q = next((r for r in metrics28 if r["model"] == model and
                  r["summary"] == summary and r["period"] == "full"), None)
        add(f"T4.1_{model}_{summary}_{key.upper()}", expected,
            float(q[key]) if q else None, 0.005, 1,
            str(phase28 or "fit_metrics.csv"))
    movement = read_csv(movement28) if movement28 else []
    for model, summary, expected in (("M2", "M2", -3.268),
                                      ("PARK", "PARK_mean", -1.388)):
        q = next((r for r in movement if r["model"] == model and
                  r["summary"] == summary), None)
        add(f"T4.2_{model}_COG", expected,
            float(q["fitted_cog_movement"]) if q else None,
            0.0005, 1, str(movement28 or "movement_metrics.csv"))
        add(f"T4.2_{model}_PEAK_MOVE", -11,
            int(q["fitted_peak_movement"]) if q else None,
            0, 1, str(movement28 or "movement_metrics.csv"))
    runtime = read_csv(runtime28) if runtime28 else []
    inventory = read_csv(inventory28) if inventory28 else []
    for model, seconds, coordinates in (("M2", 608.23, 64),
                                         ("PARK", 2153.27, 269)):
        q = next((r for r in runtime if r["model"] == model), None)
        i = next((r for r in inventory if r["model"] == model), None)
        add(f"T4.3_{model}_SECONDS", seconds,
            float(q["total_seconds"]) if q else None,
            0.005, 1, str(runtime28 or "runtime_metrics.csv"))
        add(f"T4.3_{model}_COORDINATES", coordinates,
            int(i["free_coordinates"]) if i else None,
            0, 1, str(inventory28 or "parameter_inventory.csv"))
    forecast = locate(args.forecast, "metrics.csv")
    scored = locate(args.forecast, "scored_weekly.csv")
    forecast_rows = read_csv(forecast) if forecast else []
    expected_forecast = (("2024-W01-W52", 52, 1157.5, 1132.0),
                         ("2025-W01-W14", 14, 712.7, 1470.6),
                         ("2025-W15-W52", 38, 526.6, 2629.0),
                         ("2026-W01-W44", 44, 643.3, 2187.3))
    for period, weeks, m2, naive in expected_forecast:
        q = next((r for r in forecast_rows if r["period"] == period), None)
        add(f"T5.1_{period}_WEEKS", weeks,
            int(q["weeks"]) if q else None, 0, 1,
            str(forecast or "metrics.csv"))
        for key, value in (("mae_M2", m2), ("mae_naive", naive)):
            add(f"T5.1_{period}_{key}", value,
                float(q[key]) if q else None, 0.0500000001, 1,
                str(forecast or "metrics.csv"))
    add("CH5_SCORED_WEEKS", 148,
        sum(bool(r["cases_observed"]) for r in read_csv(scored)) if scored else None,
        0, 1,
        str(scored or "scored_weekly.csv"))
    fields = list(rows[0])
    with (args.out / "verification_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    counts = {status: sum(r["status"] == status for r in rows)
              for status in ("PASS", "FAIL", "BLOCKED", "NOT_RUN")}
    (args.out / "verification_summary.md").write_text(
        "# Numerical verification\n\n" +
        "\n".join(f"- {key}: {value}" for key, value in counts.items()) +
        "\n\nA PASS requires computed output at the stated reproduction level. "
        "Input-row checks only confirm a file is present; they do not verify its numerical contents.\n")
    (args.out / "manifest.json").write_text(json.dumps({
        "python": platform.python_version(), "system": platform.platform(),
        "commands": "scripts/verify_report.py", "counts": counts}, indent=2) + "\n")
    print(counts)
    raise SystemExit(1 if counts["FAIL"] else 0)


if __name__ == "__main__":
    main()
