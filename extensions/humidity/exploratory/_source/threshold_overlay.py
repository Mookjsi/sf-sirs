"""Park et al. (2025) Figure 1A를 Python/Matplotlib으로 재현한다.

입력: 바탕화면 cases.csv, 2013~2023년 RSV 발생건수.
추가 입력: 바탕화면 temperature_smooth.csv의 스무딩된 temperature 열.
원본과 같은 계산:
  1. 53주를 52주로 합치고 지역/연도/주차별 평균을 반올림한다.
  2. 각 지역의 2013~2019년 최대 주간 발생건수로 나눈다.
  3. 상대값의 색상 상한은 2로 두며, 지역은 위도 내림차순으로 배치한다.

원본 R 코드의 그림을 Python으로 옮긴 것이므로 글꼴/여백은 조금 다를 수 있다.
추가 표시를 얹을 때는 create_heatmap()이 반환하는 ax를 사용한다.
data의 plot_x는 연도+(주차-1)/52, plot_y는 지역 행의 중심 위치이다.
기온이 12/14/16℃를 통과하면 파랑, 20/22/24℃를 통과하면 빨강으로 표시한다.
주별 스무딩값 사이를 선형 보간하며, 상승·하강 통과를 모두 표시한다.
이 표식은 실외기온 기준의 통과 시점이며 실제 냉난방 가동 기록은 아니다.

실행: python paper_rsv_heatmap.py
저장만 하기: python paper_rsv_heatmap.py --no-show
기온 표식 없이 보기: python paper_rsv_heatmap.py --without-temperature --no-show
"""

import argparse
import json
import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


CASES_PATH = Path.home() / "Desktop" / "cases.csv"
TEMPERATURE_SMOOTH_PATH = Path.home() / "Desktop" / "temperature_smooth.csv"
OUTPUT_DIR = Path(__file__).resolve().parent
FIRST_YEAR, LAST_YEAR = 2013, 2023
COLOR_MAX = 2.0
REGION_ALIASES = {"Hyōgo": "Hyogo", "Naoasaki": "Nagasaki"}
THRESHOLD_COLORS = {
    12: "#66b7ff", 14: "#168dff", 16: "#0062ff",
    20: "#ff8b80", 22: "#ff4238", 24: "#f00018",
}
SOURCE_URL = (
    "https://github.com/parksw3/perturbation/blob/"
    "475c7dcdadc12987033d332b101f025625074301/figure/figure1.R"
)


def prepare_data(cases_path=CASES_PATH):
    """원본 script_data.R 및 figure1.R와 같은 계산으로 표시용 자료를 만든다."""
    raw = pd.read_csv(cases_path)
    raw = raw.loc[raw["year"].between(FIRST_YEAR, LAST_YEAR),
                  ["prefecture", "year", "week", "cases", "lat", "long"]].copy()
    if raw.empty or not raw["week"].between(1, 53).all():
        raise ValueError("No usable 2013-2023 data, or invalid week numbers.")
    if not np.isfinite(raw[["cases", "lat", "long"]].to_numpy()).all():
        raise ValueError("Cases or geographic coordinates contain missing values.")
    if (raw["cases"] < 0).any():
        raise ValueError("Case counts must be nonnegative.")

    # 원본처럼 53주를 52주로 묶어 평균을 구한다.
    raw["week"] = raw["week"].clip(upper=52)
    data = raw.groupby(["prefecture", "year", "week"], as_index=False).agg(
        cases=("cases", "mean"), lat=("lat", "first"), long=("long", "first"))
    # R의 round와 동일한 짝수 반올림: 2.5→2, 3.5→4.
    data["cases"] = np.rint(data["cases"]).astype(np.int64)
    baseline = (data.loc[data["year"] < 2020].groupby("prefecture")["cases"]
                .max().rename("baseline_max"))
    if (baseline <= 0).any():
        raise ValueError("Each prefecture needs a positive pre-2020 maximum.")
    data = data.merge(baseline, on="prefecture", validate="many_to_one")
    data["relative_cases_uncapped"] = data["cases"] / data["baseline_max"]
    data["relative_cases"] = data["relative_cases_uncapped"].clip(upper=COLOR_MAX)

    order = (data.groupby("prefecture")["lat"].mean()
             .sort_values(ascending=False, kind="stable").index.tolist())
    expected = pd.MultiIndex.from_product(
        [order, range(FIRST_YEAR, LAST_YEAR + 1), range(1, 53)],
        names=["prefecture", "year", "week"])
    actual = pd.MultiIndex.from_frame(data[["prefecture", "year", "week"]])
    if len(order) != 47 or len(expected.difference(actual)) or len(actual.difference(expected)):
        raise ValueError("Expected complete weekly data for 47 prefectures, 2013-2023.")
    row_positions = {name: index for index, name in enumerate(order)}
    data["plot_x"] = data["year"] + (data["week"] - 1) / 52
    data["plot_y"] = data["prefecture"].map(row_positions)
    return data.sort_values(["plot_y", "year", "week"]).reset_index(drop=True)


def create_heatmap(data):
    """fig, ax, data를 반환한다. 추가할 선이나 점은 ax 위에 그릴 수 있다."""
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    # 논문 Fig.1A처럼 큰 제목 없이 A 패널만 그린다.
    fig = plt.figure(figsize=(10.8, 6.3))
    fig.canvas.manager.set_window_title("Park 2025: RSV heatmap (Figure 1A)")
    ax = fig.add_axes([0.13, 0.245, 0.855, 0.69])
    order = data.drop_duplicates("plot_y").sort_values("plot_y")["prefecture"].tolist()
    matrix = data.pivot(index="plot_y", columns="plot_x", values="relative_cases").sort_index()
    x_centers = matrix.columns.to_numpy()
    # 원본 geom_raster와 같은 셀 중심 위치와 너비를 사용한다.
    x_edges = np.r_[x_centers - 0.5 / 52, x_centers[-1] + 0.5 / 52]
    y_edges = np.arange(len(order) + 1) - 0.5
    image = ax.pcolormesh(
        x_edges, y_edges, matrix.to_numpy(), shading="flat", cmap="viridis",
        norm=Normalize(0, COLOR_MAX), rasterized=True, edgecolors="none",
    )
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlim(x_edges[0], LAST_YEAR + 1)
    ax.set_yticks(np.arange(len(order)), order)
    ax.set_xticks(np.arange(FIRST_YEAR, LAST_YEAR + 2))
    ax.set_xlabel("Year", fontsize=8, labelpad=6)
    ax.tick_params(axis="y", labelsize=7, length=1.5, width=0.45, pad=2, colors="#333333")
    ax.tick_params(axis="x", labelsize=7.5, length=2, width=0.5, pad=3, colors="#333333")
    for year in range(FIRST_YEAR, LAST_YEAR + 2):
        ax.axvline(year, color="white", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
    for spine in ax.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.5)

    colorbar_ax = fig.add_axes([0.55, 0.105, 0.12, 0.037])
    colorbar = fig.colorbar(image, cax=colorbar_ax, orientation="horizontal")
    colorbar.set_ticks([0, 1, 2], labels=["0", "1", ">2"])
    colorbar.ax.tick_params(labelsize=7, length=2, width=0.5, direction="in",
                           color="white", labelcolor="black", pad=3)
    colorbar.outline.set_visible(False)
    colorbar.ax.text(-0.10, 0.5, "Relative cases", transform=colorbar.ax.transAxes,
                     ha="right", va="center", fontsize=8, color="black")
    return fig, ax, data


def find_threshold_crossings(x, temperature, threshold):
    """주별 값 사이의 교차 위치와 양쪽 관측 인덱스를 반환한다.

    기준값에 닿았다가 같은 쪽으로 돌아오면 통과로 세지 않는다.
    정확히 기준값인 관측이 연속되면 반대편으로 이동한 경우에만 한 번 세며,
    기준값에 처음 도달한 주차를 대표 위치로 사용한다.
    """
    x = np.asarray(x, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    signed = temperature - threshold
    nonzero = np.flatnonzero(signed != 0)
    crossings = []
    for left, right in zip(nonzero[:-1], nonzero[1:]):
        if np.sign(signed[left]) == np.sign(signed[right]):
            continue
        if right == left + 1:
            fraction = (threshold - temperature[left]) / (temperature[right] - temperature[left])
            position = x[left] + fraction * (x[right] - x[left])
        else:
            position = x[left + 1]
        direction = "warming" if signed[right] > signed[left] else "cooling"
        crossings.append((float(position), int(left), int(right), direction))
    return crossings


def prepare_temperature_crossings(heatmap_data, temperature_path=TEMPERATURE_SMOOTH_PATH):
    """도도부현별 스무딩 기온에서 6개 기준온도의 통과 지점을 계산한다."""
    temperature_path = Path(temperature_path)
    if not temperature_path.exists() and temperature_path == TEMPERATURE_SMOOTH_PATH:
        temperature_path = OUTPUT_DIR.parent / "temperature_smoothing" / "temperature_smooth.csv"
    raw = pd.read_csv(temperature_path)
    required = {"key", "year", "week", "temperature"}
    if not required.issubset(raw.columns):
        raise ValueError("Smoothed temperature CSV needs key, year, week, temperature.")
    raw = raw.copy()
    raw["prefecture"] = raw["key"].replace(REGION_ALIASES)
    if raw.duplicated(["prefecture", "year", "week"]).any():
        raise ValueError("Duplicate prefecture/year/week in smoothed temperatures.")
    if not raw["week"].between(1, 52).all() or not np.isfinite(raw["temperature"]).all():
        raise ValueError("Invalid weeks or missing smoothed temperatures.")
    row_positions = heatmap_data.drop_duplicates("prefecture").set_index("prefecture")["plot_y"].to_dict()
    if set(raw["prefecture"]) != set(row_positions):
        raise ValueError("Temperature regions must match all 47 heatmap prefectures.")

    records = []
    ranges = []
    for prefecture, group in raw.groupby("prefecture"):
        ordered = group.sort_values(["year", "week"]).reset_index(drop=True)
        x = ((ordered["year"] - FIRST_YEAR) * 52 + ordered["week"] - 1).to_numpy(dtype=float)
        y = ordered["temperature"].to_numpy(dtype=float)
        if len(x) < 2 or not np.all(np.diff(x) == 1):
            raise ValueError(f"Temperature weeks must be consecutive for {prefecture}.")
        ranges.append((x[0], x[-1]))
        for threshold in THRESHOLD_COLORS:
            for position, left, right, direction in find_threshold_crossings(x, y, threshold):
                plot_x = FIRST_YEAR + position / 52
                if not (FIRST_YEAR <= plot_x < LAST_YEAR + 1):
                    continue
                year_offset = int(np.floor(position / 52))
                records.append({
                    "prefecture": prefecture, "threshold_degC": threshold,
                    "threshold_group": "heating" if threshold < 20 else "cooling",
                    "direction": direction, "crossing_year": FIRST_YEAR + year_offset,
                    "crossing_week": position - year_offset * 52 + 1,
                    "plot_x": plot_x, "plot_y": row_positions[prefecture],
                    "left_year": int(ordered.loc[left, "year"]),
                    "left_week": int(ordered.loc[left, "week"]),
                    "left_temperature_degC": float(y[left]),
                    "right_year": int(ordered.loc[right, "year"]),
                    "right_week": int(ordered.loc[right, "week"]),
                    "right_temperature_degC": float(y[right]),
                })
    columns = ["prefecture", "threshold_degC", "threshold_group", "direction", "crossing_year",
               "crossing_week", "plot_x", "plot_y", "left_year", "left_week", "left_temperature_degC",
               "right_year", "right_week", "right_temperature_degC"]
    crossings = pd.DataFrame.from_records(records, columns=columns).sort_values(
        ["plot_y", "plot_x", "threshold_degC"]).reset_index(drop=True)
    if len(set(ranges)) != 1:
        raise ValueError("Temperature coverage must be identical across prefectures.")
    start, end = ranges[0]
    start_year, end_year = FIRST_YEAR + int(start // 52), FIRST_YEAR + int(end // 52)
    start_week, end_week = int(start % 52) + 1, int(end % 52) + 1
    metadata = {
        "temperature_source": str(temperature_path.resolve()), "temperature_column": "temperature",
        "temperature_units": "degrees Celsius", "interpolation": "linear between weekly smoothed values",
        "thresholds_degC": list(THRESHOLD_COLORS), "colors": THRESHOLD_COLORS,
        "directions": ["warming", "cooling"], "number_of_crossings": len(crossings),
        "coverage_start": {"year": start_year, "week": start_week},
        "coverage_end": {"year": end_year, "week": end_week},
        "extrapolation": False, "region_aliases": REGION_ALIASES,
        "interpretation": "Outdoor temperature threshold crossings; proxy, not observed HVAC operation",
        "coverage_label": f"{start_year} W{start_week:02d} - {end_year} W{end_week:02d}",
    }
    return crossings, metadata


def add_temperature_overlay(fig, ax, crossings, metadata):
    """각 지역의 행 안에 짧은 세로 표식을 겹쳐 그린다."""
    handles = []
    for threshold, color in THRESHOLD_COLORS.items():
        selected = crossings.loc[crossings["threshold_degC"].eq(threshold)]
        ax.vlines(selected["plot_x"], selected["plot_y"] - 0.30,
                  selected["plot_y"] + 0.30, color=color, linewidth=0.6, alpha=0.5, zorder=5)
        handles.append(Line2D([], [], color=color, linestyle="none", marker="|",
                              markersize=8, markeredgewidth=1.2, alpha=0.5, label=f"{threshold}°C"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.56, 0.035), ncol=6,
               frameon=False, fontsize=7.5, columnspacing=1.5, handletextpad=0.3,
               title="Smoothed temperature crossings: blue = heating thresholds, red = cooling thresholds",
               title_fontsize=7.5)
    fig.text(0.13, 0.012,
             f"Temperature coverage: {metadata['coverage_label']} only | Upward and downward crossings | Weekly linear interpolation",
             fontsize=6.5, color="#555555")


def save_outputs(fig, data, output_dir, cases_path=CASES_PATH, crossings=None, temperature_metadata=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = "paper_rsv_heatmap_temperature_crossings" if crossings is not None else "paper_rsv_heatmap"
    for extension in ["png", "svg", "pdf"]:
        fig.savefig(output_dir / f"{stem}.{extension}", dpi=300, facecolor="white")
    data.to_csv(output_dir / "heatmap_data.csv", index=False)
    metadata = {
        "source_code": SOURCE_URL, "input_cases": str(Path(cases_path).resolve()),
        "plot": "Park et al. 2025 Figure 1A, Python reproduction",
        "normalization": "cases / prefecture maximum during 2013-2019; capped at 2",
        "week_53": "Combine with week 52 by mean, then round half to even as in R",
        "prefectures": data["prefecture"].nunique(), "weeks_per_prefecture": 572,
        "region_order_top_to_bottom": data.drop_duplicates("plot_y")["prefecture"].tolist(),
        "period": [FIRST_YEAR, LAST_YEAR], "color_limits": [0, COLOR_MAX],
        "colormap": "viridis",
    }
    if temperature_metadata is not None:
        metadata["temperature_overlay"] = temperature_metadata
    if crossings is not None:
        crossings.to_csv(output_dir / "temperature_threshold_crossings.csv", index=False, encoding="utf-8-sig")
        regions = data.drop_duplicates("plot_y")["prefecture"].tolist()
        count_index = pd.MultiIndex.from_product(
            [regions, list(THRESHOLD_COLORS), ["warming", "cooling"]],
            names=["prefecture", "threshold_degC", "direction"])
        counts = crossings.groupby(["prefecture", "threshold_degC", "direction"]).size()
        counts.reindex(count_index, fill_value=0).rename("n_crossings").reset_index().to_csv(
            output_dir / "temperature_threshold_counts.csv", index=False, encoding="utf-8-sig")
    (output_dir / "heatmap_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Saved 47 prefectures x 572 weeks to {output_dir.resolve()}")
    return output_dir / f"{stem}.png"


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
    parser.add_argument("--cases-path", type=Path, default=CASES_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--temperature-path", type=Path, default=TEMPERATURE_SMOOTH_PATH)
    parser.add_argument("--without-temperature", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    interactive = configure_display(args.no_show)
    fig, ax, data = create_heatmap(prepare_data(args.cases_path))

    crossings, temperature_metadata = None, None
    if not args.without_temperature:
        # 원래 히트맵도 별도 이름으로 저장한다.
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for extension in ["png", "svg", "pdf"]:
            fig.savefig(args.output_dir / f"paper_rsv_heatmap_baseline.{extension}", dpi=300, facecolor="white")
        crossings, temperature_metadata = prepare_temperature_crossings(data, args.temperature_path)
        add_temperature_overlay(fig, ax, crossings, temperature_metadata)
        print(f"Temperature threshold crossings: {len(crossings)}")

    figure_path = save_outputs(fig, data, args.output_dir, args.cases_path, crossings, temperature_metadata)
    if args.no_show:
        plt.close(fig)
    elif interactive:
        plt.show(block=True)
    else:
        plt.close(fig)
        os.startfile(str(figure_path.resolve()))


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
