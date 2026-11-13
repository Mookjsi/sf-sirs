"""일본 5개 지역의 대표 지점에서 연도별 주별 습도를 한 창에 표시한다.

입력: 바탕화면 humidity.csv, 비습 kg/kg를 g/kg로 변환.
자료 단위는 도시가 아닌 도도부현이다. 각 도시가 속한 도도부현을 사용한다.
홋카이도: 삿포로(홋카이도), 혼슈: 도쿄, 시코쿠: 고치,
규슈: 후쿠오카, 류큐: 나하(오키나와).
연도별 선은 오래된 연도부터 연회색→검은색이다.
2020년은 자료가 있는 1~9주까지만 그린다.

실행: python five_region_humidity.py
저장만 하기: python five_region_humidity.py --no-show
"""

import argparse
import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd


HUMIDITY_PATH = Path.home() / "Desktop" / "humidity.csv"
OUTPUT_DIR = Path(__file__).resolve().parent
REGIONS = [
    ("Hokkaido", "Hokkaido", "Sapporo (Hokkaido)"),
    ("Honshu", "Tokyo", "Tokyo"),
    ("Shikoku", "Kochi", "Kochi"),
    ("Kyushu", "Fukuoka", "Fukuoka"),
    ("Ryukyu", "Okinawa", "Naha (Okinawa)"),
]


def load_humidity():
    raw = pd.read_csv(HUMIDITY_PATH)
    prefectures = [prefecture for _, prefecture, _ in REGIONS]
    data = raw.loc[raw["key"].isin(prefectures), ["key", "year", "week", "humidity"]].copy()
    if set(data["key"]) != set(prefectures):
        raise ValueError("At least one representative prefecture is missing.")
    if data.duplicated(["key", "year", "week"]).any():
        raise ValueError("Duplicate prefecture/year/week observations.")
    if not data["week"].between(1, 52).all():
        raise ValueError("Week numbers must be between 1 and 52.")
    if not np.isfinite(data["humidity"]).all() or (data["humidity"] < 0).any():
        raise ValueError("Humidity contains missing, non-finite, or negative values.")
    data["humidity_g_kg"] = data["humidity"] * 1000
    return data.sort_values(["key", "year", "week"]).reset_index(drop=True)


def draw_panels(data, output_dir):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, axes = plt.subplots(3, 2, figsize=(14, 10.5), sharex=True, sharey=True)
    fig.canvas.manager.set_window_title("Japan: humidity in five regions")
    years = sorted(data["year"].unique())
    shades = np.linspace(0.76, 0.0, len(years))
    ymax = float(np.ceil(data["humidity_g_kg"].max() / 5) * 5)
    legend_handles = []

    for panel_index, (ax, (region, prefecture, location)) in enumerate(zip(axes.flat, REGIONS)):
        local = data.loc[data["key"].eq(prefecture)]
        for year, shade in zip(years, shades):
            weekly = local.loc[local["year"].eq(year)].set_index("week")["humidity_g_kg"]
            # 관측하지 않은 주차는 선을 이어 만들지 않도록 NaN으로 둔다.
            weekly = weekly.reindex(range(1, 53))
            line, = ax.plot(weekly.index, weekly.values, color=str(shade), linewidth=1.55,
                            label=str(year))
            if panel_index == 0:
                legend_handles.append(line)
        ax.set_title(f"{chr(65 + panel_index)}   {region}: {location}", loc="left", fontsize=12, pad=10)
        ax.set_xlim(1, 52)
        ax.set_ylim(0, ymax)
        ax.set_xticks([1, 10, 20, 30, 40, 52])
        ax.set_yticks(np.arange(0, ymax + 1, 5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.tick_params(axis="x", labelbottom=True)
        ax.tick_params(which="minor", length=2, color="#aaaaaa")
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#e7e7e7", linewidth=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#999999")

    legend_ax = axes.flat[-1]
    legend_ax.set_axis_off()
    legend_ax.legend(handles=legend_handles, title="Year", loc="upper left", ncol=4,
                     frameon=False, fontsize=10, title_fontsize=12,
                     handlelength=2.4, labelspacing=1.0, columnspacing=1.8,
                     bbox_to_anchor=(0.04, 1.0))
    legend_ax.text(0.07, 0.26, "Earlier years: light gray\nLater years: black\n2020: weeks 1-9 only",
                   transform=legend_ax.transAxes, color="#555555", fontsize=10, linespacing=1.7)
    fig.suptitle("Weekly humidity in five representative locations", x=0.075, y=0.975,
                 ha="left", fontsize=19)
    fig.text(0.075, 0.93, "Prefecture-level data | 2013-2020 | Same axes in all panels",
             color="#555555", fontsize=11)
    fig.supxlabel("Week of year", y=0.025, fontsize=12)
    fig.supylabel("Specific humidity (g/kg)", x=0.018, fontsize=12)
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.08, top=0.875, hspace=0.34, wspace=0.15)

    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "five_region_humidity.png", dpi=220, facecolor="white")
    fig.savefig(output_dir / "five_region_humidity.svg", facecolor="white")
    print(f"Saved 5 panels with {len(years)} annual curves each to {output_dir.resolve()}")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    interactive = not args.no_show
    if args.no_show:
        plt.switch_backend("Agg")
    else:
        try:
            # 그래프 창을 지원하는지 먼저 확인한다.
            probe = plt.figure()
            plt.close(probe)
        except Exception as exc:
            if type(exc).__module__ != "_tkinter":
                raise
            # Tcl/Tk가 빠진 Python에서는 저장한 그림을 기본 이미지 창으로 연다.
            plt.switch_backend("Agg")
            interactive = False
    draw_panels(load_humidity(), args.output_dir)
    if args.no_show:
        plt.close("all")
    elif interactive:
        plt.show(block=True)
    else:
        plt.close("all")
        os.startfile(str((args.output_dir / "five_region_humidity.png").resolve()))


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
