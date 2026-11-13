"""Reproduce collaborator exploratory analyses with explicit inputs and outputs."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COMPONENTS = ("variance", "tokyo", "temperature_smoothing", "temperature_heatmap",
              "humidity_heatmap", "humidity_profiles", "threshold_overlay")
EXPECTED_COLUMNS = {
    "cases": {"prefecture", "year", "week", "time", "cases", "lat", "long"},
    "temperature": {"key", "year", "week", "temperature"},
    "humidity": {"key", "year", "week", "humidity"},
    "population": {"prefecture", "pop", "island"},
}
ALIASES = {"Hyōgo": "Hyogo", "Naoasaki": "Nagasaki"}


def source(name):
    return importlib.import_module(f"extensions.humidity.exploratory._source.{name}")


def validate_inputs(input_dir):
    tables, hashes = {}, {}
    for name, columns in EXPECTED_COLUMNS.items():
        path = input_dir / f"{name}.csv"
        if not path.is_file():
            raise FileNotFoundError(f"Missing required source input: {path}")
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        table = pd.read_csv(path)
        if not columns.issubset(table) or table.empty:
            raise ValueError(f"{path.name}: missing columns or empty input")
        tables[name] = table
    cases, temperature, humidity, population = (tables[n] for n in EXPECTED_COLUMNS)
    for name in ("cases", "temperature", "humidity"):
        table = tables[name]
        key = "prefecture" if name == "cases" else "key"
        if table.duplicated([key, "year", "week"]).any() or table[[key,"year","week"]].isna().any().any():
            raise ValueError(f"{name}: duplicate or missing prefecture-week key")
        if table[key].nunique() != 47:
            raise ValueError(f"{name}: expected 47 prefectures")
        limit = 53 if name == "cases" else 52
        if not table.week.between(1, limit).all():
            raise ValueError(f"{name}: invalid calendar week")
        measure = {"cases":"cases", "temperature":"temperature", "humidity":"humidity"}[name]
        values = pd.to_numeric(table[measure], errors="raise")
        if not np.isfinite(values).all():
            raise ValueError(f"{name}: missing or nonfinite measurement")
        if name == "cases" and ((values < 0) | (values % 1 != 0)).any():
            raise ValueError("cases: expected nonnegative integer counts")
        if name == "humidity" and ((values < 0) | (values > 0.1)).any():
            raise ValueError("humidity: expected specific humidity in kg/kg")
        if name == "temperature" and ((values < -50) | (values > 60)).any():
            raise ValueError("temperature: expected degrees Celsius")
    case_names = set(cases.prefecture)
    for name in ("temperature", "humidity"):
        if {ALIASES.get(v,v) for v in tables[name].key} != case_names:
            raise ValueError(f"{name}: prefecture mapping differs from cases")
    if population.prefecture.nunique() != 47:
        raise ValueError("population: expected 47 prefectures")
    population_names = {v.rsplit("-", 1)[0] for v in population.prefecture}
    if population_names != case_names:
        raise ValueError("population: prefecture aliases differ from cases")
    if len(temperature) != 17531 or len(humidity) != 17531:
        raise ValueError("climate: expected 17,531 observations in each source")
    for name in ("temperature", "humidity"):
        climate = tables[name]
        if (climate.year.min(),climate.year.max(),climate.loc[climate.year.eq(2020),"week"].max()) != (2013,2020,9):
            raise ValueError(f"{name}: expected 2013 W01 through 2020 W09")
    if (cases.year.min(),cases.year.max()) != (2013,2023):
        raise ValueError("cases: expected 2013 through 2023")
    return tables, hashes


def save_jpg(name, out):
    plt.gcf().savefig(out / f"{name}.jpg", dpi=300, facecolor="white")
    plt.close("all")


def run_variance(tables, input_dir, out):
    mod = source("case_variance")
    data = mod.load_data(input_dir / "cases.csv")
    results, matrix, detail = mod.calculate_variances(data)
    ratios = mod.calculate_variance_ratios(matrix, results)
    mod.plot_heatmap(matrix, out, ratios)
    save_jpg("variance_heatmap", out)
    mod.make_report(results, matrix, detail, out, ratios)
    return {"source":"codes/var.py", "tables":["all_prefecture_variances.csv","variance_matrix.csv","variance_ratios.csv"]}


def run_tokyo(tables, out):
    data = tables["cases"]
    tokyo = data.loc[data.prefecture.eq("Tokyo") & data.year.between(2013,2021)].sort_values("time")
    if tokyo.empty:
        raise ValueError("Tokyo observations are missing")
    fig, ax = plt.subplots(figsize=(14,5))
    ax.plot(tokyo.time,tokyo.cases,color="steelblue",linewidth=1.5)
    starts = tokyo.groupby("year").time.min()
    ax.set_xticks(starts.values,starts.index)
    peaks = tokyo[tokyo.week.between(15,52)].loc[lambda x:x.groupby("year").cases.idxmax()]
    for row in peaks[peaks.year.ne(2020)].itertuples():
        ax.axvspan(row.time-20,row.time+20,color="dimgray",alpha=.35)
        ax.axvline(row.time,color="black",linestyle="--")
    ax.set(xlabel="Year",ylabel="Reported RSV cases per week",title="Tokyo reported RSV cases, 2013–2021")
    fig.tight_layout(); fig.savefig(out/"tokyo_peak_windows.jpg",dpi=300);plt.close(fig)
    tokyo.to_csv(out/"tokyo_observations.csv",index=False)
    return {"source":"codes/Tokyo_plot.py", "correction":"Selected Tokyo; original selected Tottori", "tables":["tokyo_observations.csv"]}


def run_smoothing(tables, input_dir, out):
    mod=source("temperature_smoothing")
    result,curves,diag=mod.smooth_temperature(input_dir/"temperature.csv")
    fig=mod.plot_comparison(curves)
    result.to_csv(out/"temperature_smooth.csv",index=False)
    diag.to_csv(out/"smoothing_diagnostics.csv",index=False)
    fig.savefig(out/"temperature_smoothing.jpg",dpi=300);plt.close(fig)
    return {"source":"codes/temperature_smoothing/temperature_smoothing.py", "tables":["temperature_smooth.csv","smoothing_diagnostics.csv"]}


def run_heatmap(kind, tables, input_dir, out):
    mod=source(kind)
    setattr(mod,"CASES_PATH",input_dir/"cases.csv")
    setattr(mod,"HUMIDITY_PATH",input_dir/"humidity.csv")
    setattr(mod,"TEMPERATURE_PATH",input_dir/"temperature.csv")
    data,matrix,first=mod.load_humidity() if kind=="humidity_heatmap" else mod.load_temperature()
    mod.draw_heatmap(data,matrix,first,out)
    save_jpg(kind,out)
    data.to_csv(out/f"{kind}_data.csv",index=False)
    return {"source":f"codes/{kind}/{kind}.py", "tables":[f"{kind}_data.csv"]}


def run_profiles(tables,input_dir,out):
    mod=source("humidity_profiles");mod.HUMIDITY_PATH=input_dir/"humidity.csv"
    data=mod.load_humidity();mod.draw_panels(data,out);save_jpg("humidity_profiles",out)
    data.to_csv(out/"humidity_profiles_data.csv",index=False)
    return {"source":"codes/five_region_humidity/five_region_humidity.py","tables":["humidity_profiles_data.csv"]}


def run_overlay(tables,input_dir,out):
    smooth=out/"temperature_smooth.csv"
    if not smooth.exists():
        raise FileNotFoundError("Threshold overlay requires temperature_smooth.csv in the run output")
    mod=source("threshold_overlay")
    data=mod.prepare_data(input_dir/"cases.csv")
    fig,ax,data=mod.create_heatmap(data)
    crossings,metadata=mod.prepare_temperature_crossings(data,smooth)
    mod.add_temperature_overlay(fig,ax,crossings,metadata)
    fig.savefig(out/"rsv_temperature_thresholds.jpg",dpi=300);plt.close(fig)
    data.to_csv(out/"rsv_heatmap_data.csv",index=False)
    crossings.to_csv(out/"temperature_threshold_crossings.csv",index=False)
    return {"source":"codes/paper_rsv_heatmap/paper_rsv_heatmap.py", "tables":["rsv_heatmap_data.csv","temperature_threshold_crossings.csv"]}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir",type=Path,default=Path("data/climate/source"))
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--analysis",choices=("all",)+COMPONENTS,default="all")
    args=p.parse_args(argv)
    if args.out.exists():
        p.error("--out must be a fresh directory")
    tables,hashes=validate_inputs(args.input_dir)
    args.out.mkdir(parents=True)
    (args.out/"input_hashes.json").write_text(json.dumps(hashes,indent=2)+"\n")
    selected=COMPONENTS if args.analysis=="all" else (args.analysis,)
    if "threshold_overlay" in selected and "temperature_smoothing" not in selected:
        selected=("temperature_smoothing",)+selected
    provenance={}
    for name in selected:
        try:
            if name=="variance": result=run_variance(tables,args.input_dir,args.out)
            elif name=="tokyo": result=run_tokyo(tables,args.out)
            elif name=="temperature_smoothing": result=run_smoothing(tables,args.input_dir,args.out)
            elif name in ("temperature_heatmap","humidity_heatmap"):result=run_heatmap(name,tables,args.input_dir,args.out)
            elif name=="humidity_profiles":result=run_profiles(tables,args.input_dir,args.out)
            else:result=run_overlay(tables,args.input_dir,args.out)
            provenance[name]={"status":"PASS",**result}
        except Exception as exc:
            provenance[name]={"status":"FAIL","error":str(exc)}
            break
    (args.out/"provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    (args.out/"environment.json").write_text(json.dumps({"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,"matplotlib":matplotlib.__version__,"scipy":importlib.import_module("scipy").__version__},indent=2)+"\n")
    (args.out/"configuration.json").write_text(json.dumps({"analysis":args.analysis,"selected":selected},indent=2)+"\n")
    (args.out/"validation.json").write_text(json.dumps({"inputs":"PASS","analyses":provenance},indent=2)+"\n")
    if any(v["status"]!="PASS" for v in provenance.values()):
        raise SystemExit(1)

if __name__=="__main__":
    main()
