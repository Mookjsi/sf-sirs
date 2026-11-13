"""Annual peak-centered RSV timing variance (population variance, 41 weeks).

Run: python var.py
Optional: python var.py --no-show --output-dir PATH
Requires: pandas, numpy, matplotlib
"""

from pathlib import Path
import argparse
import html
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


CSV_PATH = Path(r"C:\Users\손님\Desktop\cases.csv")
YEARS = list(range(2013, 2022))
HALF_WINDOW = 20
EXCLUDED_YEAR = 2020
VARIANCE_RATIO_THRESHOLD = 0.8
TOTAL_CASES_RATIO_THRESHOLD = 1.2
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "variance_outputs"


def load_data(path):
    df = pd.read_csv(path)
    required = ["prefecture", "year", "week", "time", "cases", "lat"]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if df[required].isna().any().any():
        raise ValueError("Missing values in required columns; do not treat missing cases as zero.")
    for column in ["year", "week", "time", "cases", "lat"]:
        df[column] = pd.to_numeric(df[column], errors="raise")
    if (df["cases"] < 0).any() or (df["cases"] % 1 != 0).any():
        raise ValueError("Cases must be nonnegative integer counts.")
    if (df["time"] % 1 != 0).any():
        raise ValueError("time must be an integer weekly index.")
    if df.duplicated(["prefecture", "time"]).any():
        raise ValueError("Duplicate weekly time index within a prefecture.")
    return df.sort_values(["prefecture", "time"]).reset_index(drop=True)


def calculate_variances(df):
    # 논문의 reorder(prefecture, lat)와 같은 위도순: 북쪽이 위, 남쪽이 아래.
    # 47개 도도부현을 각각 계산한다.
    order = (
        df.groupby("prefecture")["lat"].mean()
        .sort_values(ascending=False, kind="stable").index.tolist()
    )
    records = []
    tokyo_2013 = None
    for prefecture in order:
        region = df[df["prefecture"] == prefecture].sort_values("time")
        for year in YEARS:
            # 최대 케이스 주차는 해당 연도의 15~52주(양 끝 포함)에서만 찾는다.
            year_data = region[(region["year"] == year) & region["week"].between(15, 52)]
            if year_data.empty:
                raise ValueError(f"No data in weeks 15-52 for {prefecture}, {year}")
            # 동률이면 시간순으로 먼저 나온 주를 선택한다.
            peak = year_data.loc[year_data["cases"].idxmax()]
            peak_time = int(peak["time"])
            start, end = peak_time - HALF_WINDOW, peak_time + HALF_WINDOW
            # 연도 필터 전의 전체 자료에서 선택: 12월/1월 경계도 연속 41주.
            window = region[region["time"].between(start, end)].copy()
            expected = np.arange(start, end + 1)
            if not np.array_equal(window["time"].to_numpy(), expected):
                raise ValueError(f"Incomplete 41-week window: {prefecture}, {year}")

            x = np.arange(1, 2 * HALF_WINDOW + 2, dtype=float)
            counts = window["cases"].to_numpy(dtype=float)
            total = int(counts.sum())
            sum_cx = float(np.dot(counts, x))
            sum_cx2 = float(np.dot(counts, x ** 2))
            mean = sum_cx / total if total else np.nan
            weighted_ss = counts * (x - mean) ** 2
            # 각 주차 x를 cases번 반복한 표본의 모분산과 정확히 같은 계산.
            variance = float(weighted_ss.sum() / total) if total else np.nan
            first, last = window.iloc[0], window.iloc[-1]
            records.append({
                "prefecture": prefecture,
                "year": year,
                "peak_week": int(peak["week"]),
                "peak_cases": int(peak["cases"]),
                "peak_time": peak_time,
                "window_start_year": int(first["year"]),
                "window_start_week": int(first["week"]),
                "window_end_year": int(last["year"]),
                "window_end_week": int(last["week"]),
                "n_weeks": len(window),
                "total_cases": total,
                "sum_cases_x": sum_cx,
                "sum_cases_x2": sum_cx2,
                "mean_window_week": mean,
                "population_variance_weeks2": variance,
                "status": "ok" if total else "undefined_zero_cases",
            })
            if prefecture == "Tokyo" and year == 2013:
                window["window_week"] = x.astype(int)
                window["offset_from_peak"] = window["time"] - peak_time
                window["cases_x"] = counts * x
                window["cases_x2"] = counts * x ** 2
                window["mean_window_week"] = mean
                window["squared_deviation"] = (x - mean) ** 2
                window["cases_squared_deviation"] = weighted_ss
                window["variance_contribution"] = weighted_ss / total
                columns = [
                    "window_week", "offset_from_peak", "year", "week", "time", "cases",
                    "cases_x", "cases_x2", "mean_window_week", "squared_deviation",
                    "cases_squared_deviation", "variance_contribution",
                ]
                tokyo_2013 = window[columns].reset_index(drop=True)
    if tokyo_2013 is None:
        raise ValueError("Tokyo 2013 data not found")
    results = pd.DataFrame(records)
    matrix = results.pivot(
        index="prefecture", columns="year", values="population_variance_weeks2"
    ).reindex(index=order, columns=YEARS)
    return results, matrix, tokyo_2013


def calculate_variance_ratios(matrix, results):
    # 2020년을 건너뛴 순서로 비교하므로 2021년의 이전 연도는 2019년이다.
    comparison_years = [year for year in matrix.columns if year != EXCLUDED_YEAR]
    totals = results.set_index(["prefecture", "year"])["total_cases"]
    records = []
    for previous_year, year in zip(comparison_years, comparison_years[1:]):
        for prefecture in matrix.index:
            previous = float(matrix.loc[prefecture, previous_year])
            current = float(matrix.loc[prefecture, year])
            # 분모가 0이거나 분산이 결측이면 비율을 정의하지 않고 표시하지 않는다.
            ratio = current / previous if np.isfinite(current) and np.isfinite(previous) and previous > 0 else np.nan
            previous_cases = int(totals.loc[(prefecture, previous_year)])
            current_cases = int(totals.loc[(prefecture, year)])
            cases_ratio = current_cases / previous_cases if previous_cases > 0 else np.nan
            # 분산 감소와 같은 피크 구간의 총 케이스 증가를 모두 만족해야 한다.
            red_boundary = bool(
                np.isfinite(ratio) and ratio < VARIANCE_RATIO_THRESHOLD
                and np.isfinite(cases_ratio) and cases_ratio >= TOTAL_CASES_RATIO_THRESHOLD
            )
            records.append({
                "prefecture": prefecture,
                "previous_year": previous_year,
                "year": year,
                "previous_variance": previous,
                "current_variance": current,
                "variance_ratio": ratio,
                "previous_total_cases": previous_cases,
                "current_total_cases": current_cases,
                "total_cases_ratio": cases_ratio,
                "red_boundary": red_boundary,
            })
    return pd.DataFrame(records)


def plot_heatmap(matrix, output_dir, comparisons):
    # 모든 지역/연도에 동일한 선형 색상 척도. 지역별 정규화하지 않는다.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(12, 12.5))
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("#d9dde3")
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if not len(finite):
        raise ValueError("No defined variances for the heatmap")
    image = ax.imshow(
        np.ma.masked_invalid(values), aspect="auto", interpolation="nearest",
        cmap=cmap, norm=Normalize(vmin=0, vmax=max(float(finite.max()), 1.0)),
    )
    ax.set_xticks(range(len(YEARS)), YEARS)
    ax.set_yticks(range(len(matrix)), matrix.index, fontsize=8.5)
    ax.set_xlabel("Year of annual peak", labelpad=10)
    ax.set_title("RSV timing variance around each annual peak", loc="left", fontsize=16, pad=36)
    ax.text(
        0, 1.018, "47 prefectures | peak +/-20 weeks (41 weeks) | population variance",
        transform=ax.transAxes, fontsize=10, color="#59616c", va="bottom",
    )
    ax.set_xticks(np.arange(-0.5, len(YEARS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(matrix), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.45, alpha=0.3)
    # 2020년은 계산표에는 남기되 히트맵에서는 열 전체를 검은색으로 덮는다.
    if EXCLUDED_YEAR in matrix.columns:
        excluded_col = matrix.columns.get_loc(EXCLUDED_YEAR)
        ax.axvspan(excluded_col - 0.5, excluded_col + 0.5, color="black", linewidth=0, zorder=3)
    # 지역 한 행 높이의 경계선만 표시한다. 2021/2019는 2021 칸 왼쪽에 표시.
    for comparison in comparisons.itertuples(index=False):
        if comparison.red_boundary:
            row = matrix.index.get_loc(comparison.prefecture)
            col = matrix.columns.get_loc(comparison.year)
            ax.vlines(col - 0.5, row - 0.5, row + 0.5, color="#ff2020", linewidth=2.5, zorder=5)
    ax.tick_params(which="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax, orientation="horizontal", fraction=0.035, pad=0.065, aspect=38)
    colorbar.set_label("Population variance (weeks²)   |   yellow: smaller / purple: larger", labelpad=8)
    colorbar.outline.set_visible(False)
    fig.legend(
        handles=[
            Patch(facecolor="black", label="2020 excluded"),
            Line2D([0], [0], color="#ff2020", linewidth=2.5,
                   label=f"Variance ratio < {VARIANCE_RATIO_THRESHOLD} and total cases ratio >= {TOTAL_CASES_RATIO_THRESHOLD}"),
        ],
        loc="lower center", bbox_to_anchor=(0.56, 0.012), ncol=2,
        frameon=False, fontsize=9,
    )
    fig.subplots_adjust(left=0.155, right=0.97, top=0.925, bottom=0.085)
    fig.savefig(output_dir / "variance_heatmap.png", dpi=240, facecolor="white")
    fig.savefig(output_dir / "variance_heatmap.svg", facecolor="white")
    return fig


def make_report(results, matrix, detail, output_dir, comparisons):
    tokyo = results[results["prefecture"] == "Tokyo"].copy()
    tokyo["window"] = (
        tokyo["window_start_year"].astype(str) + "-W" + tokyo["window_start_week"].astype(str).str.zfill(2)
        + " ~ " + tokyo["window_end_year"].astype(str) + "-W" + tokyo["window_end_week"].astype(str).str.zfill(2)
    )
    tokyo_table = tokyo[[
        "year", "peak_week", "peak_cases", "window", "total_cases", "mean_window_week",
        "population_variance_weeks2",
    ]]
    row = tokyo[tokyo["year"] == 2013].iloc[0]
    n, sx, sx2 = int(row["total_cases"]), row["sum_cases_x"], row["sum_cases_x2"]
    mean, variance = row["mean_window_week"], row["population_variance_weeks2"]
    ss = detail["cases_squared_deviation"].sum()
    explanation = (
        "Tokyo 2013: population variance of case occurrence week\n"
        f"Peak: 2013-W{int(row['peak_week']):02d}, cases={int(row['peak_cases'])}\n"
        f"Window: {row['window']} (41 consecutive weeks, includes adjacent-year observations)\n"
        "x = window week 1..41; peak is x=21. c_x = reported cases in that week.\n"
        "Equivalent expanded sample: repeat x exactly c_x times.\n"
        f"N = sum(c_x) = {n}\n"
        f"sum(c_x*x) = {sx:.0f}\n"
        f"sum(c_x*x^2) = {sx2:.0f}\n"
        f"mean = sum(c_x*x)/N = {sx:.0f}/{n} = {mean:.10f}\n"
        f"sum(c_x*(x-mean)^2) = {ss:.10f}\n"
        f"population variance = {ss:.10f}/{n} = {variance:.10f} weeks^2\n"
        f"Cross-check: {sx2:.0f}/{n} - ({sx:.0f}/{n})^2 = {sx2/n-(sx/n)**2:.10f}\n"
        "The variance is unchanged if x is shifted to -20..20 or the global time index.\n"
        "This is variance of occurrence timing, not variance of the 41 case-count values.\n"
    )
    results.to_csv(output_dir / "all_prefecture_variances.csv", index=False, encoding="utf-8-sig")
    matrix.to_csv(output_dir / "variance_matrix.csv", encoding="utf-8-sig")
    comparisons.to_csv(output_dir / "variance_ratios.csv", index=False, encoding="utf-8-sig")
    tokyo_table.to_csv(output_dir / "tokyo_variances.csv", index=False, encoding="utf-8-sig")
    detail.to_csv(output_dir / "tokyo_2013_calculation.csv", index=False, encoding="utf-8-sig")
    (output_dir / "tokyo_2013_calculation.txt").write_text(
        explanation + "\n41-week calculation table:\n" + detail.to_string(index=False, float_format=lambda x: f"{x:.8f}"),
        encoding="utf-8",
    )
    display_tokyo = tokyo_table.rename(columns={
        "year": "연도", "peak_week": "최대 주차", "peak_cases": "최대 케이스",
        "window": "41주 구간", "total_cases": "총 케이스 N", "mean_window_week": "평균 주차 (1~41)",
        "population_variance_weeks2": "모분산 (주²)",
    })
    display_detail = detail.rename(columns={
        "window_week": "구간 주차 x", "offset_from_peak": "피크 대비 주차", "year": "연도",
        "week": "실제 주차", "cases": "케이스 c", "cases_x": "c×x", "cases_x2": "c×x²",
        "mean_window_week": "평균 μ", "squared_deviation": "(x−μ)²",
        "cases_squared_deviation": "c×(x−μ)²", "variance_contribution": "분산 기여도",
    })
    page = """<!doctype html><html lang="ko"><meta charset="utf-8"><title>RSV 주차 모분산</title>
<style>body{font:15px/1.65 'Malgun Gothic',sans-serif;color:#233044;max-width:1200px;margin:36px auto;padding:0 22px}
h1,h2{line-height:1.3}table{border-collapse:collapse;font-size:13px;width:100%;font-variant-numeric:tabular-nums}
th,td{padding:8px 10px;border-bottom:1px solid #e1e6ed;text-align:right;white-space:nowrap}th{background:#edf2f7}
tr:nth-child(even){background:#f7f9fb}.scroll{overflow:auto}pre{background:#f0f4f8;padding:18px;white-space:pre-wrap}
img{width:100%;max-width:950px}a{color:#1765a0}</style>
<h1>RSV 발생 주차의 모분산</h1>
<p>2013~2021년(2020년 포함), 47개 도도부현. 매년 15~52주 중 최대 케이스 주차를 중심으로 연속 41주를 사용합니다.
연도 경계를 넘는 구간도 실제 관측자료를 이어서 사용합니다. 최대값 동률은 먼저 나온 주를 선택합니다.</p>
<p>케이스 수 c를 주차 x의 빈도로 사용합니다. μ = Σcx / Σc, 모분산 = Σc(x−μ)² / Σc.
구간 주차는 1~41이며 피크는 21입니다. 작은 분산은 케이스가 특정 시점에 더 모였음을 뜻합니다.</p>
<h2>지역별·연도별 히트맵</h2><p>논문 Figure 1A처럼 북쪽에서 남쪽으로 정렬합니다.
전체 칸에 동일한 선형 색상 척도를 사용하며 노란색은 작은 분산, 보라색은 큰 분산입니다.
2020년 열은 모든 지역에서 검은색으로 표시합니다. 각 지역의 분산을 이전 비교 연도 분산으로 나눈 값이
0.8 미만이고, 같은 피크 ±20주 구간의 총 케이스 수를 이전 비교 연도의 구간 내 총 케이스 수로 나눈 값이
1.2 이상일 때만 해당 연도 칸의 왼쪽 경계를 빨간색으로 표시합니다. 비교에서 2020년을 건너뛰므로
2021년은 2019년과 비교하고 빨간선은 2021년 칸 왼쪽에 표시합니다. 계산표에는 2020년 값을 남깁니다.</p>
<img src="variance_heatmap.png" alt="47개 지역 9개 연도의 주차 모분산 히트맵">
<h2>도쿄: 연도별 결과</h2>"""
    page += '<div class="scroll">' + display_tokyo.to_html(index=False, border=0, float_format=lambda x: f"{x:.6f}") + '</div>'
    display_ratios = comparisons[comparisons["prefecture"] == "Tokyo"].drop(columns="prefecture").rename(columns={
        "previous_year": "이전 비교 연도", "year": "연도",
        "previous_variance": "이전 분산", "current_variance": "현재 분산",
        "variance_ratio": "현재 / 이전 분산", "red_boundary": "빨간 경계선",
        "previous_total_cases": "이전 총 케이스", "current_total_cases": "현재 총 케이스",
        "total_cases_ratio": "현재 / 이전 총 케이스",
    })
    page += '<h2>도쿄: 분산 비율과 경계선</h2><div class="scroll">' + display_ratios.to_html(index=False, border=0, float_format=lambda x: f"{x:.6f}") + '</div>'
    page += '<h2>도쿄 2013년: 계산 과정</h2><pre>' + html.escape(explanation) + '</pre>'
    page += '<div class="scroll">' + display_detail.to_html(index=False, border=0, float_format=lambda x: f"{x:.6f}") + '</div>'
    page += '<p>각 행의 분산 기여도를 합하면 위의 모분산이 됩니다. N=0인 구간은 분산을 정의할 수 없어 결측으로 표시합니다.</p></html>'
    (output_dir / "variance_report.html").write_text(page, encoding="utf-8")
    print("\nTokyo annual results:\n")
    print(tokyo_table.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    print("\n" + explanation)
    print(detail.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    print(f"\nSaved results: {output_dir.resolve()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=CSV_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    results, matrix, detail = calculate_variances(load_data(args.csv))
    comparisons = calculate_variance_ratios(matrix, results)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_heatmap(matrix, args.output_dir, comparisons)
    make_report(results, matrix, detail, args.output_dir, comparisons)
    if args.no_show:
        plt.close("all")
    else:
        plt.show()


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
