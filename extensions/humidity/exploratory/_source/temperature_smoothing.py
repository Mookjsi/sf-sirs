"""일본 주별 기온 자료를 스무딩하고 5개 대표 지역에서 비교한다.

입력: 바탕화면 temperature.csv, 섭씨 기온.
방법: 각 지역의 전체 시계열에 3차 스무딩 스플라인을 적합한다.
      평활 강도는 GCV로 자동 선택한다(lam=None).
출력: temperature_smooth.csv의 temperature는 스무딩값,
      temperature_original은 원자료이다. 원본 파일은 수정하지 않는다.
      원래 행과 주차를 유지하고, 곡선을 관측 범위 밖으로 외삽하지 않는다.

이 코드는 과거 시계열의 흐름을 매끈하게 요약하는 회고적 분석이다.
실행: python temperature_smoothing.py
저장만 하기: python temperature_smoothing.py --no-show

방법 문서:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_smoothing_spline.html
"""

import argparse
import json
import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import make_smoothing_spline
import numpy as np
import pandas as pd


TEMPERATURE_PATH = Path.home() / "Desktop" / "temperature.csv"
OUTPUT_DIR = Path(__file__).resolve().parent
REPRESENTATIVES = [
    ("Hokkaido", "Hokkaido: Sapporo (Hokkaido)"),
    ("Tokyo", "Honshu: Tokyo"),
    ("Kochi", "Shikoku: Kochi"),
    ("Fukuoka", "Kyushu: Fukuoka"),
    ("Okinawa", "Ryukyu: Naha (Okinawa)"),
]


def smooth_temperature(input_path=TEMPERATURE_PATH):
    """원래 행 순서를 유지한 스무딩 자료, 지역별 곡선, 진단표를 반환한다."""
    raw = pd.read_csv(input_path)
    required = {"key", "year", "week", "temperature"}
    if not required.issubset(raw.columns) or raw.empty:
        raise ValueError("temperature.csv must contain key, year, week, temperature.")
    if raw.duplicated(["key", "year", "week"]).any():
        raise ValueError("Duplicate prefecture/year/week observations.")
    if raw[["key", "year", "week"]].isna().any().any():
        raise ValueError("Region or time identifiers contain missing values.")
    if not raw["week"].between(1, 52).all():
        raise ValueError("This dataset must use weeks 1-52.")
    if not np.isfinite(raw["temperature"]).all():
        raise ValueError("Temperature contains missing or non-finite values.")

    result = raw.copy()
    result["temperature_original"] = raw["temperature"]
    curves, diagnostics = {}, []
    first_year = int(raw["year"].min())
    for region, group in raw.groupby("key", sort=True):
        ordered = group.sort_values(["year", "week"])
        x = ((ordered["year"] - first_year) * 52 + ordered["week"] - 1).to_numpy(dtype=float)
        y = ordered["temperature"].to_numpy(dtype=float)
        if len(x) < 5 or not np.all(np.diff(x) == 1):
            raise ValueError(f"Expected a continuous weekly series for {region}.")
        # 연말에 곡선을 끊지 않고 해당 지역의 모든 주차를 한 번에 적합한다.
        spline = make_smoothing_spline(x, y, lam=None)
        fitted = spline(x)
        if not np.isfinite(fitted).all():
            raise ValueError(f"Non-finite smoothed temperatures for {region}.")
        result.loc[ordered.index, "temperature"] = fitted
        curves[region] = {"x": x, "original": y, "smooth": fitted,
                          "spline": spline, "first_year": first_year}
        diagnostics.append({
            "key": region, "n_weeks": len(x),
            "rmse_change_degC": float(np.sqrt(np.mean((fitted - y) ** 2))),
            "max_abs_change_degC": float(np.max(np.abs(fitted - y))),
            "mean_original_degC": float(y.mean()),
            "mean_smooth_degC": float(fitted.mean()),
            "min_original_degC": float(y.min()), "max_original_degC": float(y.max()),
            "min_smooth_degC": float(fitted.min()), "max_smooth_degC": float(fitted.max()),
            "roughness_original": float(np.mean(np.diff(y, n=2) ** 2)),
            "roughness_smooth": float(np.mean(np.diff(fitted, n=2) ** 2)),
        })
    if not {name for name, _ in REPRESENTATIVES}.issubset(curves):
        raise ValueError("Some representative prefectures are missing.")
    return result, curves, pd.DataFrame(diagnostics)


def plot_comparison(curves):
    """원자료와 연속 스무딩 곡선을 한 창의 5개 세로 패널에 겹쳐 그린다."""
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(5, 1, figsize=(15, 11), sharex=True, sharey=True)
    fig.canvas.manager.set_window_title("Temperature: original and smoothed, five regions")
    first_year = curves[REPRESENTATIVES[0][0]]["first_year"]
    xmin = min(c["x"][0] for c in curves.values())
    xmax = max(c["x"][-1] for c in curves.values())
    values = np.concatenate([np.r_[c["original"], c["smooth"]] for c in curves.values()])
    ymin = float(np.floor(values.min() / 5) * 5)
    ymax = float(np.ceil(values.max() / 5) * 5)

    for i, (ax, (region, title)) in enumerate(zip(axes, REPRESENTATIVES)):
        curve = curves[region]
        x = curve["x"]
        plot_x = first_year + x / 52
        # 원래 주별 관측은 회색 점과 선으로, 스무딩은 빨간 곡선으로 구분한다.
        ax.plot(plot_x, curve["original"], color="#606060", linewidth=0.75,
                marker="o", markersize=1.6, markeredgewidth=0, zorder=1)
        # 저장한 주별 값과 같은 스플라인을 관측 범위 안에서만 촘촘하게 평가한다.
        dense_x = np.linspace(x[0], x[-1], int((x[-1] - x[0]) * 12) + 1)
        ax.plot(first_year + dense_x / 52, curve["spline"](dense_x),
                color="#cf3038", linewidth=1.155, zorder=2)
        ax.set_title(f"{chr(65 + i)}   {title}", loc="left", fontsize=10.5, pad=5)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(np.arange(ymin, ymax + 1, 10))
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#e7e7e7", linewidth=0.6)
        for year in range(first_year + 1, int(first_year + xmax / 52) + 1):
            ax.axvline(year, color="#dddddd", linewidth=0.6, linestyle=(0, (3, 3)), zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#999999")
        ax.tick_params(labelsize=9)

    axes[-1].set_xlim(first_year + xmin / 52, first_year + xmax / 52)
    axes[-1].set_xticks(np.arange(first_year, int(first_year + xmax / 52) + 1))
    axes[-1].set_xlabel("Year", labelpad=8)
    fig.supylabel("Weekly mean temperature (°C)", x=0.025, fontsize=12)
    fig.suptitle("Weekly temperature: original and smoothed", x=0.075, y=0.98,
                 ha="left", fontsize=19)
    fig.text(0.075, 0.944, "Prefecture-level data | 2013 W01 - 2020 W09 | Cubic smoothing spline (GCV)",
             fontsize=10, color="#555555")
    handles = [
        Line2D([], [], color="#606060", linewidth=1, marker="o", markersize=3, label="Original weekly data"),
        Line2D([], [], color="#cf3038", linewidth=1.26, label="Smoothed curve"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.984),
               frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.065, top=0.90, hspace=0.34)
    return fig


def save_outputs(result, diagnostics, fig, output_dir, input_path):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "temperature_smooth.csv", index=False, encoding="utf-8-sig")
    diagnostics.to_csv(output_dir / "smoothing_diagnostics.csv", index=False, encoding="utf-8-sig")
    for extension in ["png", "svg", "pdf"]:
        fig.savefig(output_dir / f"temperature_smoothing_comparison.{extension}", dpi=240, facecolor="white")
    metadata = {
        "input": str(Path(input_path).resolve()), "units": "degrees Celsius",
        "method": "scipy.interpolate.make_smoothing_spline", "spline_degree": 3,
        "smoothing_strength": "Automatically selected separately for each prefecture using GCV (lam=None)",
        "fit_period": "Entire weekly series per prefecture; no reset at year boundaries",
        "rows": len(result), "n_prefectures": result["key"].nunique(),
        "row_order": "Same as input", "temperature_column": "Smoothed weekly values",
        "temperature_original_column": "Original weekly observations",
        "extrapolation": False,
        "analysis_type": "Retrospective smoothing; not a forecast",
        "method_documentation": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_smoothing_spline.html",
        "representative_keys": [key for key, _ in REPRESENTATIVES],
    }
    (output_dir / "smoothing_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Saved {len(result)} rows across {result['key'].nunique()} prefectures to {output_dir.resolve()}")


def configure_display(no_show):
    if no_show:
        plt.switch_backend("Agg")
        return False
    try:
        probe = plt.figure()
        plt.close(probe)
        return True
    except Exception as exc:
        if type(exc).__module__ != "_tkinter":
            raise
        plt.switch_backend("Agg")
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=TEMPERATURE_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    interactive = configure_display(args.no_show)
    result, curves, diagnostics = smooth_temperature(args.input)
    fig = plot_comparison(curves)
    save_outputs(result, diagnostics, fig, args.output_dir, args.input)
    if args.no_show:
        plt.close(fig)
    elif interactive:
        plt.show(block=True)
    else:
        plt.close(fig)
        os.startfile(str((args.output_dir / "temperature_smoothing_comparison.png").resolve()))


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
