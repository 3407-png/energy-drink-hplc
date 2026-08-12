"""외부 검량선 · 표준물 첨가법 · 변환 상수 계산 (탐구 3).

계산 흐름
    1) 증류수 바탕 검량선  A = m*C + b        -> 순수 용매에서의 응답계수 m
    2) 음료 바탕 표준물 첨가 A = k*C_add + A0  -> x절편 |−A0/k| = 바이알 내 원래 농도
    3) 희석배수(10배)를 곱해 음료 원액 농도, 1회 제공량을 곱해 mg/캔
    4) 변환 상수 f = (역산한 참 농도) / (외부 검량선으로 읽은 겉보기 농도)

읽기 전에 알아둘 것 (연구 가설 둘째와 직접 관련):
    표준물 첨가법이 보정해 주는 것은 **곱셈형(multiplicative) 매트릭스 효과**,
    즉 매트릭스 때문에 단위 농도당 응답(기울기)이 달라지는 현상이다.
    반대로 **덧셈형(additive) 간섭**, 즉 리보플라빈처럼 같은 시간대에 같이 나와서
    피크 면적을 통째로 더해 버리는 공용리(co-elution)는 표준물 첨가법으로
    보정되지 않는다. 오히려 A0 를 부풀려 x절편을 과대평가하게 만든다.
    -> 덧셈형 간섭은 수학이 아니라 크로마토그래피(분리능)로 해결해야 한다.
    이 모듈은 그 크기를 `additive_sensitivity()` 로 정량해 준다.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .config import (
    CALIB_GROUP,
    COMPOUND_ORDER,
    COMPOUNDS,
    DRINKS,
    PREP,
    QC,
    SAMPLE_GROUP,
    THRESHOLDS,
    Compound,
    Drink,
)
from .stats import LinearFit, linear_fit, propagate_ratio


# ---------------------------------------------------------------------------
# 검량선
# ---------------------------------------------------------------------------

@dataclass
class CalibrationResult:
    compound: Compound
    fit: LinearFit
    rsd_by_level: dict[float, float] = field(default_factory=dict)

    @property
    def lod_ppm(self) -> float:
        return self.fit.lod

    @property
    def loq_ppm(self) -> float:
        return self.fit.loq

    @property
    def passes_linearity(self) -> bool:
        return bool(self.fit.r2 >= QC.min_r2)


def build_calibrations(df: pd.DataFrame) -> dict[str, CalibrationResult]:
    """증류수 바탕 검량선(Vial 1~4)을 성분별로 만든다."""
    from .stats import rsd_percent

    calib = df[df["group"] == CALIB_GROUP]
    if calib.empty:
        raise ValueError(
            "검량선 데이터(group=calib)가 없습니다. "
            "외부 검량선 없이는 변환 상수를 계산할 수 없습니다."
        )

    results: dict[str, CalibrationResult] = {}
    for cmp_key in COMPOUND_ORDER:
        sub = calib[calib["compound"] == cmp_key]
        if sub.empty:
            continue
        level_mean = sub.groupby("spike_ppm")["peak_area"].mean()
        fit = linear_fit(level_mean.index.values, level_mean.values)
        rsds = {
            float(lvl): rsd_percent(g["peak_area"].values)
            for lvl, g in sub.groupby("spike_ppm")
        }
        results[cmp_key] = CalibrationResult(COMPOUNDS[cmp_key], fit, rsds)
    return results


# ---------------------------------------------------------------------------
# 표준물 첨가법
# ---------------------------------------------------------------------------

@dataclass
class StandardAdditionResult:
    drink: Drink
    compound: Compound
    fit: LinearFit

    # 외부 검량선(같은 성분)
    calibration: CalibrationResult | None = None

    # 반복 주입 정밀도
    rsd_by_level: dict[float, float] = field(default_factory=dict)

    # ---- 바이알 내 농도 (희석된 상태) ---------------------------------
    @property
    def c_vial_ppm(self) -> float:
        """|x절편|. 표준물 첨가법으로 역산한 바이알 내 참 농도."""
        return abs(self.fit.x_intercept)

    @property
    def se_c_vial_ppm(self) -> float:
        return self.fit.se_x_intercept

    def ci_c_vial_ppm(self, level: float = QC.confidence_level) -> tuple[float, float]:
        lo, hi = self.fit.x_intercept_ci(level)
        return (abs(hi), abs(lo)) if self.fit.x_intercept < 0 else (lo, hi)

    # ---- 음료 원액 농도 --------------------------------------------------
    @property
    def c_drink_ppm(self) -> float:
        """음료 원액 농도 (mg/L). 바이알 농도 x 희석배수."""
        return self.c_vial_ppm * PREP.dilution_factor

    @property
    def se_c_drink_ppm(self) -> float:
        return self.se_c_vial_ppm * PREP.dilution_factor

    @property
    def mg_per_serving(self) -> float:
        """1회 제공량(캔 하나)당 mg."""
        return self.c_drink_ppm * self.drink.serving_mL / 1000.0

    @property
    def se_mg_per_serving(self) -> float:
        return self.se_c_drink_ppm * self.drink.serving_mL / 1000.0

    # ---- 외부 검량선으로 읽은 겉보기 농도 --------------------------------
    @property
    def a0_area(self) -> float:
        """무첨가 시료(Vial A1)의 신호값 = 회귀직선의 y절편."""
        return self.fit.intercept

    @property
    def c_apparent_vial_ppm(self) -> float:
        """외부 검량선법으로 구한 겉보기 농도 (간섭 미보정)."""
        if self.calibration is None:
            return float("nan")
        return self.calibration.fit.concentration_from_area(self.a0_area)

    @property
    def c_apparent_drink_ppm(self) -> float:
        return self.c_apparent_vial_ppm * PREP.dilution_factor

    # ---- 변환 상수 -------------------------------------------------------
    @property
    def conversion_constant(self) -> float:
        """f = 참 농도 / 겉보기 농도.

        f > 1 : 외부 검량선이 실제보다 과소평가 (매트릭스가 응답을 억제)
        f < 1 : 외부 검량선이 실제보다 과대평가 (공용리 등이 면적을 부풀림)
        """
        app = self.c_apparent_vial_ppm
        if not math.isfinite(app) or app == 0:
            return float("nan")
        return self.c_vial_ppm / app

    @property
    def se_conversion_constant(self) -> float:
        if self.calibration is None:
            return float("nan")
        app = self.c_apparent_vial_ppm
        # 겉보기 농도의 오차는 검량선 기울기 오차에서 주로 온다.
        se_app = (
            abs(app) * self.calibration.fit.se_slope / abs(self.calibration.fit.slope)
            if self.calibration.fit.slope
            else float("nan")
        )
        _, se = propagate_ratio(self.c_vial_ppm, self.se_c_vial_ppm, app, se_app)
        return se

    @property
    def bias_percent(self) -> float:
        """외부 검량선법의 상대오차(%) = (겉보기 - 참) / 참 * 100."""
        if self.c_vial_ppm == 0:
            return float("nan")
        return (self.c_apparent_vial_ppm - self.c_vial_ppm) / self.c_vial_ppm * 100.0

    # ---- 매트릭스 효과 분해 ---------------------------------------------
    @property
    def slope_ratio(self) -> float:
        """k(첨가법 기울기) / m(검량선 기울기).

        1 보다 작으면 매트릭스가 응답을 억제한 것(이온화/흡광 방해),
        1 보다 크면 증강. 이것이 표준물 첨가법이 실제로 보정해 주는 성분이다.
        """
        if self.calibration is None or self.calibration.fit.slope == 0:
            return float("nan")
        return self.fit.slope / self.calibration.fit.slope

    @property
    def matrix_effect_percent(self) -> float:
        """곱셈형 매트릭스 효과 (%). 0 이면 매트릭스 영향 없음."""
        r = self.slope_ratio
        return (r - 1.0) * 100.0 if math.isfinite(r) else float("nan")

    def additive_sensitivity(
        self, fractions: tuple[float, ...] = (0.05, 0.10, 0.20)
    ) -> dict[float, float]:
        """덧셈형 공용리 간섭에 대한 민감도 분석.

        A0 중 일부가 사실 타깃이 아닌 공용리 성분(리보플라빈 등)의 면적이라면,
        x절편은 그 비율만큼 그대로 과대평가된다.
        반환값: {A0 중 간섭 기여율 -> 그때의 보정 농도(ppm)}
        """
        out = {}
        for frac in fractions:
            if self.fit.slope == 0:
                out[frac] = float("nan")
                continue
            out[frac] = abs(-(self.a0_area * (1 - frac)) / self.fit.slope)
        return out

    # ---- 품질 판정 -------------------------------------------------------
    @property
    def x_intercept_rsd_percent(self) -> float:
        if self.c_vial_ppm == 0:
            return float("nan")
        return self.se_c_vial_ppm / self.c_vial_ppm * 100.0

    def warnings(self) -> list[str]:
        w: list[str] = []
        name = f"{self.drink.name_ko} / {self.compound.name_ko}"

        if self.fit.r2 < QC.min_r2:
            w.append(
                f"[{name}] 첨가법 직선성 미달 (R2={self.fit.r2:.4f} < {QC.min_r2}). "
                "회귀 신뢰도가 낮아 x절편 역산값을 그대로 쓰기 어렵습니다."
            )
        if (
            math.isfinite(self.x_intercept_rsd_percent)
            and self.x_intercept_rsd_percent > QC.max_x_intercept_rsd_percent
        ):
            w.append(
                f"[{name}] x절편 상대표준오차 {self.x_intercept_rsd_percent:.1f}% "
                f"(기준 {QC.max_x_intercept_rsd_percent}%). 외삽 거리가 너무 멉니다 "
                "- 첨가 농도 범위를 시료 농도에 맞춰 조정하세요."
            )
        for lvl, r in self.rsd_by_level.items():
            if math.isfinite(r) and r > QC.max_rsd_percent:
                w.append(
                    f"[{name}] +{lvl:.0f} ppm 바이알 반복주입 RSD {r:.1f}% "
                    f"(기준 {QC.max_rsd_percent}%). 주입/여과 재현성을 확인하세요."
                )
        if self.calibration is not None and self.c_vial_ppm < self.calibration.loq_ppm:
            w.append(
                f"[{name}] 역산 농도 {self.c_vial_ppm:.2f} ppm 가 정량한계"
                f"({self.calibration.loq_ppm:.2f} ppm) 미만입니다. "
                "정량값으로 보고하지 말고 '검출' 수준으로만 서술하세요."
            )
        if math.isfinite(self.matrix_effect_percent) and abs(self.matrix_effect_percent) > 10:
            w.append(
                f"[{name}] 곱셈형 매트릭스 효과 {self.matrix_effect_percent:+.1f}% "
                "- 외부 검량선법을 쓰면 안 되는 수준입니다(첨가법 채택 근거)."
            )
        return w


def run_standard_addition(
    df: pd.DataFrame, calibrations: dict[str, CalibrationResult]
) -> list[StandardAdditionResult]:
    """음료 3종 x 성분 3종에 대해 표준물 첨가법 회귀를 수행한다."""
    from .stats import rsd_percent

    samples = df[df["group"] == SAMPLE_GROUP]
    if samples.empty:
        raise ValueError("음료 시료 데이터(group=sample)가 없습니다.")

    results: list[StandardAdditionResult] = []
    for drink_key in DRINKS:
        sub_d = samples[samples["sample"] == drink_key]
        if sub_d.empty:
            continue
        for cmp_key in COMPOUND_ORDER:
            sub = sub_d[sub_d["compound"] == cmp_key]
            if sub.empty:
                continue
            level_mean = sub.groupby("spike_ppm")["peak_area"].mean()
            if len(level_mean) < 3:
                continue
            fit = linear_fit(level_mean.index.values, level_mean.values)
            rsds = {
                float(lvl): rsd_percent(g["peak_area"].values)
                for lvl, g in sub.groupby("spike_ppm")
            }
            results.append(
                StandardAdditionResult(
                    drink=DRINKS[drink_key],
                    compound=COMPOUNDS[cmp_key],
                    fit=fit,
                    calibration=calibrations.get(cmp_key),
                    rsd_by_level=rsds,
                )
            )
    return results


# ---------------------------------------------------------------------------
# 라벨 대조 및 안전성 평가 (연구 목적 셋째)
# ---------------------------------------------------------------------------

@dataclass
class LabelComparison:
    drink: Drink
    measured_mg: float
    se_mg: float
    label_mg: float | None

    @property
    def recovery_percent(self) -> float:
        if not self.label_mg:
            return float("nan")
        return self.measured_mg / self.label_mg * 100.0

    @property
    def difference_mg(self) -> float:
        if self.label_mg is None:
            return float("nan")
        return self.measured_mg - self.label_mg

    @property
    def within_tolerance(self) -> bool:
        """식품 표시 허용 오차 관례상 표시량의 80~120% 이내면 일치로 본다."""
        r = self.recovery_percent
        return bool(math.isfinite(r) and 80.0 <= r <= 120.0)


def compare_with_label(results: list[StandardAdditionResult]) -> list[LabelComparison]:
    out = []
    for r in results:
        if r.compound.key != "caffeine":
            continue
        out.append(
            LabelComparison(
                drink=r.drink,
                measured_mg=r.mg_per_serving,
                se_mg=r.se_mg_per_serving,
                label_mg=r.drink.label_caffeine_mg,
            )
        )
    return out


def safety_assessment(results: list[StandardAdditionResult]) -> list[dict]:
    """청소년 기준 섭취 위해도 평가 (결론 작성용)."""
    rows = []
    for r in results:
        row = {
            "drink": r.drink.name_ko,
            "compound": r.compound.name_ko,
            "conc_ppm": r.c_drink_ppm,
            "mg_per_serving": r.mg_per_serving,
        }
        if r.compound.key == "caffeine":
            limit = THRESHOLDS.teen_daily_caffeine_mg
            row["daily_limit_mg"] = limit
            row["percent_of_limit"] = r.mg_per_serving / limit * 100.0
            row["cans_to_reach_limit"] = (
                limit / r.mg_per_serving if r.mg_per_serving > 0 else float("inf")
            )
            row["high_caffeine_label"] = r.c_drink_ppm >= THRESHOLDS.high_caffeine_ppm
        elif r.compound.adi_mg_per_kg is not None:
            limit = r.compound.adi_mg_per_kg * THRESHOLDS.teen_bw_kg
            row["daily_limit_mg"] = limit
            row["percent_of_limit"] = r.mg_per_serving / limit * 100.0
            row["cans_to_reach_limit"] = (
                limit / r.mg_per_serving if r.mg_per_serving > 0 else float("inf")
            )
            row["high_caffeine_label"] = None
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# 실험 설계 점검 (실험 전에 돌려볼 것)
# ---------------------------------------------------------------------------

def check_spike_design(expected_vial_ppm: dict[str, float]) -> list[str]:
    """첨가 농도 범위가 시료 농도에 대해 적절한지 진단한다.

    표준물 첨가법의 통용 지침: 최고 첨가 농도가 시료 농도의 약 1~3배가 되도록
    잡는다. 시료 농도에 비해 첨가량이 너무 작으면 외삽 거리가 짧아 보이지만
    기울기 정보가 부족하고, 너무 크면 원래 신호가 묻혀 x절편 오차가 커진다.
    """
    msgs: list[str] = []
    spikes = [s for s in PREP.spike_levels_ppm if s > 0]
    smax = max(spikes)

    for cmp_key, c_vial in expected_vial_ppm.items():
        c = COMPOUNDS[cmp_key]
        if c_vial <= 0:
            continue
        ratio = smax / c_vial
        line = (
            f"[{c.name_ko}] 희석 후 예상 농도 {c_vial:.1f} ppm, "
            f"최고 첨가 {smax:.0f} ppm (배율 {ratio:.1f}배)"
        )
        if ratio < 0.5:
            line += "  -> 첨가량이 너무 적습니다. 첨가 농도를 올리세요."
        elif ratio > 5:
            line += (
                "  -> 첨가량이 과도합니다. 시료 고유 신호가 묻혀 x절편 오차가 "
                "커집니다. 이 성분만 첨가 농도를 낮추거나 희석배수를 줄이세요."
            )
        else:
            line += "  -> 적정 범위."
        msgs.append(line)
    return msgs
