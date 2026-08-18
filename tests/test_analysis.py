"""분석 수식이 실제로 맞는지 검증.

핵심 검증은 '정답을 아는 가짜 데이터를 넣고 그 정답이 나오는가'다.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from hplc import analysis as an
from hplc import dataio
from hplc.config import COMPOUNDS, DRINKS, PREP
from hplc.simulate import (
    SimulationSettings,
    TRUE_DRINK_PPM,
    simulate_peak_areas,
    truth_table,
)
from hplc.stats import linear_fit, resolution, rsd_percent


# ---------------------------------------------------------------------------
# 회귀 기본
# ---------------------------------------------------------------------------

def test_linear_fit_recovers_exact_line():
    x = np.array([0.0, 20.0, 40.0, 60.0])
    y = 3.0 * x + 12.0
    fit = linear_fit(x, y)
    assert fit.slope == pytest.approx(3.0)
    assert fit.intercept == pytest.approx(12.0)
    assert fit.r2 == pytest.approx(1.0)
    assert fit.x_intercept == pytest.approx(-4.0)


def test_x_intercept_is_negative_native_concentration():
    """표준물 첨가법: 시료에 25 ppm 이 원래 있었다면 x절편은 -25."""
    native = 25.0
    k = 1500.0
    x = np.array([0.0, 20.0, 40.0, 60.0])
    y = k * (native + x)
    fit = linear_fit(x, y)
    assert abs(fit.x_intercept) == pytest.approx(native)


def test_fit_requires_three_points():
    with pytest.raises(ValueError, match="최소 3점"):
        linear_fit([0.0, 20.0], [1.0, 2.0])


def test_fit_rejects_constant_x():
    with pytest.raises(ValueError, match="동일"):
        linear_fit([10.0, 10.0, 10.0], [1.0, 2.0, 3.0])


def test_x_intercept_se_grows_with_noise():
    rng = np.random.default_rng(0)
    x = np.array([0.0, 20.0, 40.0, 60.0])
    base = 1500.0 * (25.0 + x)
    quiet = linear_fit(x, base + rng.normal(0, 100, 4))
    noisy = linear_fit(x, base + rng.normal(0, 8000, 4))
    assert noisy.se_x_intercept > quiet.se_x_intercept


def test_rsd_and_resolution():
    assert rsd_percent([100.0, 100.0, 100.0]) == pytest.approx(0.0)
    assert rsd_percent([90.0, 110.0]) == pytest.approx(14.142, rel=1e-3)
    # 두 피크 RT 3.5/5.5, 폭 0.2/0.3 -> Rs = 2*2.0/0.5 = 8
    assert resolution(3.5, 0.2, 5.5, 0.3) == pytest.approx(8.0)


def test_lod_is_smaller_than_loq():
    x = np.array([0.0, 20.0, 40.0, 60.0])
    y = 1000.0 * x + np.array([50.0, -30.0, 20.0, -10.0])
    fit = linear_fit(x, y)
    assert 0 < fit.lod < fit.loq


# ---------------------------------------------------------------------------
# 시뮬레이션 왕복 검증
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def analyzed():
    df = simulate_peak_areas(SimulationSettings(seed=1234))
    cals = an.build_calibrations(df)
    results = an.run_standard_addition(df, cals)
    return cals, results


def test_all_drink_compound_pairs_present(analyzed):
    _, results = analyzed
    # 시료/성분 구성이 바뀌어도 따라오도록 config 에서 기대값을 만든다
    assert len(results) == len(DRINKS) * len(COMPOUNDS)


def test_calibration_linearity(analyzed):
    cals, _ = analyzed
    for cal in cals.values():
        # 4점 검량선 + 주입 재현성 1.2% 조건에서 현실적으로 기대할 수 있는 수준
        assert cal.fit.r2 > 0.998
        assert cal.passes_linearity


def test_standard_addition_recovers_truth_without_additive_interference():
    """덧셈형 간섭이 없으면 x절편 역산은 참값을 거의 정확히 되찾아야 한다."""
    st = SimulationSettings(seed=7, include_additive=False)
    df = simulate_peak_areas(st)
    cals = an.build_calibrations(df)
    results = an.run_standard_addition(df, cals)

    for r in results:
        true_ppm = TRUE_DRINK_PPM[r.drink.key][r.compound.key]
        recovery = r.c_drink_ppm / true_ppm * 100
        assert 95.0 <= recovery <= 105.0, (
            f"{r.drink.key}/{r.compound.key} 회수율 {recovery:.1f}%"
        )


def test_standard_addition_beats_external_calibration():
    """곱셈형 매트릭스 효과가 있을 때, 첨가법이 외부 검량선보다 정확해야 한다.

    매트릭스 효과가 5% 이상인 조합에서는 개별적으로도 첨가법이 이겨야 하고,
    전체 평균 오차도 첨가법 쪽이 작아야 한다.
    """
    from hplc.simulate import MATRIX_RESPONSE

    st = SimulationSettings(seed=99, include_additive=False, include_matrix=True)
    df = simulate_peak_areas(st)
    cals = an.build_calibrations(df)
    results = an.run_standard_addition(df, cals)

    err_sa, err_ext = [], []
    for r in results:
        true_ppm = TRUE_DRINK_PPM[r.drink.key][r.compound.key]
        e_sa = abs(r.c_drink_ppm - true_ppm)
        e_ext = abs(r.c_apparent_drink_ppm - true_ppm)
        err_sa.append(e_sa)
        err_ext.append(e_ext)

        effect = abs(MATRIX_RESPONSE[r.drink.key][r.compound.key] - 1.0)
        if effect >= 0.05:
            assert e_sa < e_ext, f"{r.drink.key}/{r.compound.key}"

    assert np.mean(err_sa) < np.mean(err_ext)


def test_additive_interference_inflates_x_intercept():
    """덧셈형(공용리) 간섭은 표준물 첨가법으로 보정되지 않고 과대평가를 낳는다.

    연구 가설 둘째('간섭이 Y절편 상승에만 국한되므로 첨가법이 해결한다')가
    성립하지 않는 경우를 코드로 못박아 둔다. 부풀어 오르는 양은 이론적으로
    b_off / k (간섭 면적 ÷ 첨가법 기울기) 와 같아야 한다.
    """
    from hplc.simulate import ADDITIVE_OFFSET

    clean = _analyze(SimulationSettings(seed=5, include_additive=False))
    dirty = _analyze(SimulationSettings(seed=5, include_additive=True))
    clean_map = {(r.drink.key, r.compound.key): r for r in clean}

    for r in dirty:
        c = clean_map[(r.drink.key, r.compound.key)]
        assert r.c_vial_ppm > c.c_vial_ppm, f"{r.drink.key}/{r.compound.key}"

        b_off = ADDITIVE_OFFSET[r.drink.key][r.compound.key]
        predicted = b_off / r.fit.slope
        observed = r.c_vial_ppm - c.c_vial_ppm
        assert observed == pytest.approx(predicted, rel=0.15), (
            f"{r.drink.key}/{r.compound.key}: 예상 +{predicted:.3f} ppm, "
            f"실제 +{observed:.3f} ppm"
        )


def _analyze(settings):
    df = simulate_peak_areas(settings)
    cals = an.build_calibrations(df)
    return an.run_standard_addition(df, cals)


def test_slope_ratio_tracks_matrix_response():
    """기울기비(k/m)가 시뮬레이션에 넣은 매트릭스 응답계수를 되찾는가."""
    from hplc.simulate import MATRIX_RESPONSE

    results = _analyze(SimulationSettings(seed=42))
    for r in results:
        expected = MATRIX_RESPONSE[r.drink.key][r.compound.key]
        # 4점 회귀 두 개의 기울기 비이므로 ~1% 수준의 회귀 오차가 겹친다.
        assert r.slope_ratio == pytest.approx(expected, abs=0.06)


def test_conversion_constant_direction(analyzed):
    """변환 상수는 유한하고 양수여야 하며, 겉보기 농도와 곱하면 참 농도가 된다."""
    _, results = analyzed
    for r in results:
        f = r.conversion_constant
        assert math.isfinite(f) and f > 0
        assert r.c_apparent_vial_ppm * f == pytest.approx(r.c_vial_ppm, rel=1e-9)


def test_dilution_and_serving_math(analyzed):
    _, results = analyzed
    for r in results:
        assert r.c_drink_ppm == pytest.approx(r.c_vial_ppm * PREP.dilution_factor)
        assert r.mg_per_serving == pytest.approx(
            r.c_drink_ppm * r.drink.serving_mL / 1000.0
        )


def test_additive_sensitivity_reduces_estimate(analyzed):
    _, results = analyzed
    r = results[0]
    sens = r.additive_sensitivity((0.10,))
    assert sens[0.10] == pytest.approx(r.c_vial_ppm * 0.90, rel=1e-6)


# ---------------------------------------------------------------------------
# 입출력 검증
# ---------------------------------------------------------------------------

def test_template_roundtrip_is_rejected_when_empty(tmp_path):
    p = dataio.write_template(tmp_path / "t.csv")
    with pytest.raises(dataio.DataValidationError, match="채워진 행이 하나도"):
        dataio.load_peak_areas(p)


def test_simulated_csv_is_flagged_and_loads(tmp_path):
    from hplc.simulate import write_simulated_csv

    p = write_simulated_csv(tmp_path / "sim.csv")
    assert dataio.is_simulated(p)
    df = dataio.load_peak_areas(p)
    assert len(df) > 0
    assert set(df["group"]) == {"calib", "sample"}


def test_unknown_compound_is_rejected(tmp_path):
    from hplc.simulate import write_simulated_csv

    p = write_simulated_csv(tmp_path / "sim.csv")
    text = p.read_text(encoding="utf-8-sig").replace("caffeine", "kafein")
    p.write_text(text, encoding="utf-8-sig")
    with pytest.raises(dataio.DataValidationError, match="알 수 없는 값"):
        dataio.load_peak_areas(p)


def test_negative_area_is_rejected(tmp_path):
    from hplc.simulate import write_simulated_csv

    p = write_simulated_csv(tmp_path / "sim.csv")
    lines = p.read_text(encoding="utf-8-sig").splitlines()
    body = [l for l in lines if not l.startswith("#")]
    body[1] = body[1].rsplit(",", 1)[0] + ",-500"
    p.write_text("\n".join([l for l in lines if l.startswith("#")] + body),
                 encoding="utf-8-sig")
    with pytest.raises(dataio.DataValidationError, match="음수"):
        dataio.load_peak_areas(p)


def test_missing_calibration_raises():
    df = simulate_peak_areas(SimulationSettings(seed=3))
    only_samples = df[df["group"] == "sample"]
    with pytest.raises(ValueError, match="외부 검량선이 없습니다"):
        an.build_calibrations(only_samples)


def test_manual_calibration_fills_in_when_no_calib_rows(monkeypatch):
    """검량선 행이 없어도 config 에 계수를 넣어 두면 변환 상수까지 계산된다."""
    from hplc.simulate import RESPONSE_FACTOR

    monkeypatch.setattr(
        an, "EXTERNAL_CALIBRATION",
        {k: (v, 0.0) for k, v in RESPONSE_FACTOR.items()},
    )
    df = simulate_peak_areas(SimulationSettings(seed=11))
    only_samples = df[df["group"] == "sample"]

    cals = an.build_calibrations(only_samples)
    assert all(c.source == "manual" for c in cals.values())
    assert all(c.fit is None for c in cals.values())
    assert all(not math.isfinite(c.r2) for c in cals.values())

    results = an.run_standard_addition(only_samples, cals)
    assert len(results) == len(DRINKS) * len(COMPOUNDS)
    for r in results:
        assert math.isfinite(r.conversion_constant)
        assert math.isfinite(r.c_drink_ppm)


def test_truth_table_shape():
    t = truth_table()
    assert len(t) == len(DRINKS) * len(COMPOUNDS)
    assert {"true_drink_ppm", "matrix_response", "additive_offset_area"} <= set(t.columns)
