"""주별 습도 히트맵: 기온 히트맵과 같은 지역·시간 배치.

입력: 바탕화면 humidity.csv, cases.csv(지역별 위도 확인용)
습도는 비습으로, 원본 kg/kg에 1000을 곱해 g/kg로 표시한다.
색상: 건조할수록 회색, 습할수록 파랑. 모든 지역에 동일한 척도.
실행: python humidity_heatmap.py
저장만 하기: python humidity_heatmap.py --no-show
"""

from pathlib import Path
import argparse
import json
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize


DESKTOP = Path.home() / "Desktop"
HUMIDITY_PATH = DESKTOP / "humidity.csv"
CASES_PATH = DESKTOP / "cases.csv"
OUTPUT_DIR = Path(__file__).resolve().parent
REGION_ALIASES = {"Hyōgo": "Hyogo", "Naoasaki": "Nagasaki"}


def load_humidity():
    raw = pd.read_csv(HUMIDITY_PATH)
    cases = pd.read_csv(CASES_PATH)
    data = raw[["key", "year", "week", "humidity"]].copy()
    # 원본 CSV는 유지하며, 표시용 단위만 kg/kg에서 g/kg로 변환한다.
    data["humidity_g_kg"] = data["humidity"] * 1000
    data["prefecture"] = data["key"].replace(REGION_ALIASES)
    if data.duplicated(["prefecture", "year", "week"]).any():
        raise ValueError("Duplicate prefecture/year/week after region name matching.")
    if not data["week"].between(1, 52).all():
        raise ValueError("This humidity dataset uses weeks 1-52 per year.")
    if not np.isfinite(data["humidity"]).all() or (data["humidity"] < 0).any():
        raise ValueError("Humidity contains missing, non-finite, or negative values.")

    latitude = cases.groupby("prefecture")["lat"].mean()
    names = set(data["prefecture"])
    if names != set(latitude.index) or len(names) != 47:
        raise ValueError("Humidity regions and the 47 prefectures in cases.csv must match.")
    # 논문과 같은 위도순: 북쪽이 위, 남쪽이 아래.
    order = latitude.sort_values(ascending=False, kind="stable").index.tolist()
    first_year = int(data["year"].min())
    data["week_index"] = (data["year"] - first_year) * 52 + data["week"] - 1
    last_index = int(data["week_index"].max())
    matrix = data.pivot(index="prefecture", columns="week_index", values="humidity_g_kg")
    matrix = matrix.reindex(index=order, columns=range(last_index + 1))
    if matrix.isna().any().any():
        raise ValueError("A prefecture/week observation is missing; no interpolation is performed.")
    return data, matrix, first_year


def draw_heatmap(data, matrix, first_year, output_dir):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(15, 9))
    values = matrix.to_numpy(dtype=float)
    vmin = float(np.floor(values.min() / 5) * 5)
    vmax = float(np.ceil(values.max() / 5) * 5)
    if vmin == vmax:
        vmax = vmin + 5
    colors = ["#d0d0d0", "#a6c8de", "#4b93c3", "#08519c", "#08306b"]
    cmap = LinearSegmentedColormap.from_list("gray_to_blue", colors, N=256)
    cmap.set_bad("#d9dde3")
    # 실제 주별 습도값을 유지하고 연속적인 색상 척도로 표시한다.
    x_edges = first_year + np.arange(matrix.shape[1] + 1) / 52
    y_edges = np.arange(matrix.shape[0] + 1)
    image = ax.pcolormesh(
        x_edges, y_edges, values,
        cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
        shading="flat", edgecolors="none", rasterized=True,
    )
    ax.set_ylim(len(matrix), 0)
    ax.set_yticks(np.arange(len(matrix)) + 0.5, matrix.index, fontsize=8)
    last = data.sort_values(["year", "week"]).iloc[-1]
    years = np.arange(first_year, int(last["year"]) + 1)
    ax.set_xticks(years)
    ax.set_xlim(x_edges[0], x_edges[-1])
    for year in years[1:]:
        ax.axvline(year, color="white", linestyle=(0, (3, 3)), linewidth=0.8, alpha=0.7)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=3)
    ax.set_xlabel("Year", labelpad=8)
    ax.set_title("Weekly mean specific humidity across Japan", loc="left", fontsize=17, pad=28)
    ax.text(
        0, 1.015,
        f"47 prefectures, north to south | {first_year} W01 - {int(last['year'])} W{int(last['week']):02d}",
        transform=ax.transAxes, va="bottom", color="#525e6c", fontsize=10,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax, orientation="horizontal", shrink=0.53,
                            fraction=0.04, pad=0.095, aspect=42)
    colorbar.set_ticks(np.arange(vmin, vmax + 1, 5))
    colorbar.set_label("Weekly mean specific humidity (g/kg)", labelpad=7)
    colorbar.outline.set_visible(False)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.91, bottom=0.12)

    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "humidity_heatmap.png", dpi=240, facecolor="white")
    fig.savefig(output_dir / "humidity_heatmap.svg", facecolor="white")
    metadata = {
        "humidity_source": str(HUMIDITY_PATH), "latitude_source": str(CASES_PATH),
        "source_units": "kg/kg", "display_units": "g/kg", "unit_multiplier": 1000,
        "n_prefectures": len(matrix), "n_weeks": matrix.shape[1],
        "first_year": first_year, "last_year": int(last["year"]), "last_week": int(last["week"]),
        "observed_min_humidity_g_kg": float(values.min()), "observed_max_humidity_g_kg": float(values.max()),
        "color_min": vmin, "color_max": vmax, "colormap": "gray_to_blue", "color_stops": colors,
        "normalization": "shared linear scale",
        "region_aliases": REGION_ALIASES, "region_order_top_to_bottom": matrix.index.tolist(),
    }
    (output_dir / "humidity_heatmap_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Saved {len(matrix)} prefectures x {matrix.shape[1]} weeks to {output_dir.resolve()}")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    data, matrix, first_year = load_humidity()
    draw_heatmap(data, matrix, first_year, args.output_dir)
    if args.no_show:
        plt.close("all")
    else:
        plt.show(block=True)


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
