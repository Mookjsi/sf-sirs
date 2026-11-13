"""Focused tests for the independent exploratory calculations."""
from pathlib import Path
import pandas as pd
import pytest
from extensions.humidity.exploratory.runner import validate_inputs, run_tokyo
from extensions.humidity.exploratory._source import case_variance as variance
from extensions.humidity.exploratory._source import threshold_overlay as thresholds


def synthetic_tokyo():
    rows=[]
    for year in range(2013,2023):
        for week in range(1,53):
            rows.append(("Tokyo",year,week,(year-2013)*52+week,10 if week in (29,31) else 0,35.7))
    return pd.DataFrame(rows,columns=["prefecture","year","week","time","cases","lat"])


def test_peak_window_tie_variance_and_cross_year():
    results,matrix,detail=variance.calculate_variances(synthetic_tokyo())
    row=results[results.year.eq(2013)].iloc[0]
    assert (row.peak_week,row.n_weeks)==(29,41)
    assert row.population_variance_weeks2==pytest.approx(1.0)
    assert len(detail)==41
    ratios=variance.calculate_variance_ratios(matrix,results)
    assert not ratios.year.eq(2020).any()
    assert ratios[ratios.year.eq(2021)].previous_year.eq(2019).all()


def test_variance_boundary_strictness():
    matrix=pd.DataFrame({2013:[10.0],2014:[8.0]},index=["Tokyo"])
    data=pd.DataFrame([("Tokyo",2013,100),("Tokyo",2014,120)],columns=["prefecture","year","total_cases"])
    assert not variance.calculate_variance_ratios(matrix,data).red_boundary.iloc[0]
    matrix.loc["Tokyo",2014]=7.9
    assert variance.calculate_variance_ratios(matrix,data).red_boundary.iloc[0]


def test_threshold_crossing_both_directions():
    crossing=thresholds.find_threshold_crossings([0,1,2],[10,14,10],12)
    assert [x[0] for x in crossing]==[0.5,1.5]
    assert [x[3] for x in crossing]==["warming","cooling"]
    assert set(thresholds.THRESHOLD_COLORS)=={12,14,16,20,22,24}


def test_tokyo_output_uses_tokyo(tmp_path):
    run_tokyo({"cases":synthetic_tokyo()},tmp_path)
    saved=pd.read_csv(tmp_path/"tokyo_observations.csv")
    assert set(saved.prefecture)=={"Tokyo"}
    assert (tmp_path/"tokyo_peak_windows.jpg").is_file()


def test_input_validation_and_coverage():
    tables,hashes=validate_inputs(Path('data/climate/source'))
    assert len(tables['temperature'])==len(tables['humidity'])==17531
    assert tables['temperature'].year.max()==2020
    assert tables['temperature'].loc[tables['temperature'].year.eq(2020),'week'].max()==9
    assert tables['cases'].year.max()==2023
    assert len(hashes)==4
