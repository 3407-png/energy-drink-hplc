"""결과 표 · 마크다운 보고서 생성 (연구 결과 7장 초안).

여기서 나오는 표는 그대로 보고서에 옮겨 붙일 수 있도록 만들었다.
다만 모의 데이터로 돌린 경우 모든 출력물 맨 위에 경고가 박히며,
그 경고를 지우고 실측값처럼 제출하는 것은 데이터 조작이다.
"""

from __future__ import annotations

import math
from datetime import datetime
from pathlib import Path

import pandas as pd

from .analysis import (
    CalibrationResult,
    LabelComparison,
    StandardAdditionResult,
    compare_with_label,
    safety_assessment,
)
from .config import COMPOUND_ORDER, DRINKS, HPLC, PREP, QC, THRESHOLDS

SIM_WARNING = """> ## ⚠ 모의(SIMULATED) 데이터로 생성된 문서입니다
>
> 아래 수치는 HPLC로 측정한 값이 아니라 가정한 응답계수·농도로 계산한
> **예상 결과**입니다. 보고서 7장(연구 결과)에 실측값으로 옮겨 적으면
> 데이터 조작에 해당합니다. 실험 후 실측 CSV로 다시 돌려서 이 경고가
> 사라진 문서를 사용하십시오.
"""


def _f(x: float, nd: int = 2) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "-"
    return f"{x:,.{nd}f}"


def _table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(r) + " |")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# 개별 표
# ---------------------------------------------------------------------------

def calibration_table(cals: dict[str, CalibrationResult]) -> str:
    rows = []
    for key in COMPOUND_ORDER:
        cal = cals.get(key)
        if cal is None:
            continue
        r2_txt = f"{cal.r2:.5f}" if math.isfinite(cal.r2) else "-"
        if cal.source == "manual":
            verdict = "계수 직접 입력 (직선성 미평가)"
        elif cal.passes_linearity:
            verdict = "적합"
        else:
            verdict = f"미달(<{QC.min_r2})"
        rows.append([
            cal.compound.name_ko,
            f"{cal.slope:,.0f}",
            f"{cal.intercept:,.0f}",
            r2_txt,
            _f(cal.lod_ppm, 2),
            _f(cal.loq_ppm, 2),
            verdict,
        ])
    return _table(
        ["성분", "기울기 m (면적/ppm)", "y절편 b", "R²", "LOD (ppm)", "LOQ (ppm)", "직선성"],
        rows,
    )


def standard_addition_table(results: list[StandardAdditionResult]) -> str:
    rows = []
    for r in results:
        lo, hi = r.ci_c_vial_ppm()
        rows.append([
            r.drink.name_ko,
            r.compound.name_ko,
            f"{r.fit.slope:,.0f}",
            f"{r.fit.intercept:,.0f}",
            f"{r.fit.r2:.5f}",
            _f(r.fit.x_intercept, 2),
            f"{_f(r.c_vial_ppm, 2)} ± {_f(r.se_c_vial_ppm, 2)}",
            f"[{_f(lo, 2)}, {_f(hi, 2)}]",
        ])
    return _table(
        ["시료", "성분", "기울기 k", "y절편 A₀", "R²", "x절편 (ppm)",
         "바이알 농도 ± SE (ppm)", f"{int(QC.confidence_level*100)}% 신뢰구간"],
        rows,
    )


def concentration_table(results: list[StandardAdditionResult]) -> str:
    rows = []
    for r in results:
        rows.append([
            r.drink.name_ko,
            r.compound.name_ko,
            _f(r.c_vial_ppm, 2),
            _f(r.c_drink_ppm, 1),
            f"{_f(r.mg_per_serving, 1)} ± {_f(r.se_mg_per_serving, 1)}",
            f"{r.drink.serving_mL:.0f}",
        ])
    return _table(
        ["시료", "성분", "바이알 농도 (ppm)",
         f"음료 원액 농도 (ppm, ×{PREP.dilution_factor:.0f} 희석보정)",
         "1회 제공량당 함량 (mg)", "제공량 (mL)"],
        rows,
    )


def conversion_table(results: list[StandardAdditionResult]) -> str:
    rows = []
    for r in results:
        rows.append([
            r.drink.name_ko,
            r.compound.name_ko,
            _f(r.c_apparent_vial_ppm, 2),
            _f(r.c_vial_ppm, 2),
            f"{_f(r.conversion_constant, 4)} ± {_f(r.se_conversion_constant, 4)}",
            f"{r.bias_percent:+.1f}%" if math.isfinite(r.bias_percent) else "-",
            f"{r.matrix_effect_percent:+.1f}%" if math.isfinite(r.matrix_effect_percent) else "-",
        ])
    return _table(
        ["시료", "성분", "겉보기 농도 (외부 검량선)", "참 농도 (첨가법 역산)",
         "변환 상수 f", "외부 검량선 상대오차", "기울기비로 본 매트릭스 효과"],
        rows,
    )


def label_table(comparisons: list[LabelComparison]) -> str:
    rows = []
    for c in comparisons:
        rows.append([
            c.drink.name_ko,
            _f(c.label_mg, 1) if c.label_mg else "미확인",
            f"{_f(c.measured_mg, 1)} ± {_f(c.se_mg, 1)}",
            _f(c.difference_mg, 1),
            f"{c.recovery_percent:.1f}%" if math.isfinite(c.recovery_percent) else "-",
            ("일치 (80~120%)" if c.within_tolerance else "불일치")
            if c.label_mg else "라벨값 미입력",
            "" if c.drink.label_verified else "⚠ 라벨 미확인",
        ])
    return _table(
        ["시료", "라벨 표시량 (mg/캔)", "실측 역산값 (mg/캔)", "차이 (mg)",
         "표시량 대비", "판정", "비고"],
        rows,
    )


def safety_table(results: list[StandardAdditionResult]) -> str:
    rows = []
    for row in safety_assessment(results):
        limit = row.get("daily_limit_mg")
        rows.append([
            row["drink"],
            row["compound"],
            _f(row["conc_ppm"], 1),
            _f(row["mg_per_serving"], 1),
            _f(limit, 1) if limit else "-",
            f"{row['percent_of_limit']:.1f}%" if limit else "-",
            f"{row['cans_to_reach_limit']:.1f}" if limit else "-",
        ])
    return _table(
        ["시료", "성분", "농도 (ppm)", "1캔당 (mg)",
         f"청소년({THRESHOLDS.teen_bw_kg:.0f} kg) 1일 한도 (mg)",
         "1캔 = 한도의", "한도 도달 캔 수"],
        rows,
    )


def additive_sensitivity_table(results: list[StandardAdditionResult]) -> str:
    fractions = (0.05, 0.10, 0.20)
    rows = []
    for r in results:
        sens = r.additive_sensitivity(fractions)
        rows.append(
            [r.drink.name_ko, r.compound.name_ko, _f(r.c_vial_ppm, 2)]
            + [_f(sens[f], 2) for f in fractions]
        )
    return _table(
        ["시료", "성분", "보정 전 역산값 (ppm)"]
        + [f"공용리 기여 {int(f*100)}% 가정" for f in fractions],
        rows,
    )


def recovery_vs_truth_table(results: list[StandardAdditionResult]) -> str:
    """모의 데이터일 때만: 넣어 둔 정답을 얼마나 되찾았는지."""
    from .simulate import truth_table

    truth = truth_table().set_index(["drink", "compound"])
    rows = []
    for r in results:
        try:
            t = truth.loc[(r.drink.key, r.compound.key)]
        except KeyError:
            continue
        true_ppm = float(t["true_drink_ppm"])
        rec = r.c_drink_ppm / true_ppm * 100.0 if true_ppm else float("nan")
        rows.append([
            r.drink.name_ko,
            r.compound.name_ko,
            _f(true_ppm, 1),
            _f(r.c_drink_ppm, 1),
            f"{rec:.1f}%",
            _f(r.c_apparent_drink_ppm, 1),
            f"{r.c_apparent_drink_ppm / true_ppm * 100:.1f}%" if true_ppm else "-",
        ])
    return _table(
        ["시료", "성분", "입력한 참값 (ppm)", "첨가법 역산 (ppm)", "첨가법 회수율",
         "외부 검량선 (ppm)", "외부 검량선 회수율"],
        rows,
    )


# ---------------------------------------------------------------------------
# 전체 보고서
# ---------------------------------------------------------------------------

def build_report(
    cals: dict[str, CalibrationResult],
    results: list[StandardAdditionResult],
    *,
    simulated: bool,
    source: str,
    figures: list[Path] | None = None,
    extra_warnings: list[str] | None = None,
) -> str:
    parts: list[str] = []

    parts.append("# 제로 에너지 드링크 3종 HPLC 동시분석 결과")
    parts.append("")
    if simulated:
        parts.append(SIM_WARNING)
        parts.append("")
    parts.append(f"- 생성 시각: {datetime.now():%Y-%m-%d %H:%M}")
    parts.append(f"- 데이터 출처: `{source}`" + ("  **(모의 데이터)**" if simulated else ""))
    parts.append(f"- HPLC 조건: {HPLC.summary()}, 컬럼 {HPLC.column}")
    parts.append(
        f"- 전처리: 음료 원액 {PREP.drink_aliquot_mL:.0f} mL → "
        f"{PREP.final_volume_mL:.0f} mL 정용 (희석배수 {PREP.dilution_factor:.0f}), "
        f"{PREP.filter_um} µm 시린지 필터"
    )
    parts.append(
        f"- 표준물 첨가 수준: {', '.join(f'{s:.0f}' for s in PREP.spike_levels_ppm)} ppm, "
        f"시료 바이알 주입 {PREP.injections_sample}회"
    )
    parts.append("")

    parts.append("## 1. 외부 검량선 (증류수 바탕)")
    parts.append("")
    parts.append(calibration_table(cals))
    parts.append("")
    parts.append(
        "LOD/LOQ는 ICH Q2(R2)의 검량선 기반 산정식(3.3σ/S, 10σ/S)을 따랐으며, "
        "σ로는 y절편의 표준오차를 사용하였다."
    )
    parts.append("")

    parts.append("## 2. 표준물 첨가법 회귀 및 x절편 역산")
    parts.append("")
    parts.append(standard_addition_table(results))
    parts.append("")
    parts.append(
        "x절편의 표준오차는 외삽 오차식 "
        "s(x₀) = (s_y/x ÷ |k|)·√(1/n + ȳ²/(k²·Σ(xᵢ−x̄)²)) 로 계산하였다. "
        "표준물 첨가법은 데이터 구간 **바깥**을 외삽하므로, 기울기와 R²만 보고하고 "
        "이 불확도를 빼놓으면 정량값의 신뢰도를 평가할 수 없다."
    )
    parts.append("")

    parts.append("## 3. 성분별 실제 함량")
    parts.append("")
    parts.append(concentration_table(results))
    parts.append("")

    parts.append("## 4. 변환 상수 및 매트릭스 간섭 정량화")
    parts.append("")
    parts.append(conversion_table(results))
    parts.append("")
    parts.append(
        "**변환 상수 f = (첨가법으로 역산한 참 농도) / (외부 검량선으로 읽은 겉보기 농도)**\n\n"
        "- f > 1 : 매트릭스가 응답을 억제하여 외부 검량선법이 과소평가\n"
        "- f < 1 : 공용리 등으로 면적이 부풀어 외부 검량선법이 과대평가\n\n"
        "표의 마지막 열(기울기비)은 첨가법 기울기 k를 검량선 기울기 m으로 나눈 값이다. "
        "f ≈ 1/(k/m) 관계가 성립하므로, 변환 상수의 실체는 대부분 **곱셈형 매트릭스 효과의 역수**다."
    )
    parts.append("")

    parts.append("## 5. 라벨 표시량 대조 (카페인)")
    parts.append("")
    parts.append(label_table(compare_with_label(results)))
    parts.append("")

    parts.append("## 6. 섭취 위해도 평가")
    parts.append("")
    parts.append(safety_table(results))
    parts.append("")
    parts.append(
        f"청소년 카페인 권고 상한은 체중 1 kg당 {THRESHOLDS.caffeine_limit_mg_per_kg} mg "
        f"({THRESHOLDS.teen_bw_kg:.0f} kg 기준 {THRESHOLDS.teen_daily_caffeine_mg:.0f} mg/일), "
        f"고카페인 함유 표시 대상 기준은 {THRESHOLDS.high_caffeine_ppm:.0f} ppm(0.15 mg/mL)이다. "
        "아세설팜칼륨·소듐벤조에이트는 JECFA ADI를 적용하였다."
    )
    parts.append("")

    parts.append("## 7. 덧셈형(공용리) 간섭 민감도 분석")
    parts.append("")
    parts.append(
        "표준물 첨가법이 보정해 주는 것은 곱셈형 매트릭스 효과(기울기 변화)뿐이다. "
        "리보플라빈처럼 타깃과 **같은 시간대에 겹쳐 나오는** 성분이 있으면 그 면적이 "
        "A₀에 그대로 더해지고, x절편은 그만큼 과대평가된다. "
        "아래 표는 A₀ 중 일부가 타깃이 아니었다고 가정했을 때 역산값이 어떻게 변하는지 보여준다."
    )
    parts.append("")
    parts.append(additive_sensitivity_table(results))
    parts.append("")
    parts.append(
        "→ 이 간섭은 수학이 아니라 크로마토그래피로 해결해야 한다. "
        "① 피크 순도 확인(PDA 스펙트럼 비교 또는 2개 파장 면적비), "
        "② 이동상 조성·pH 조정으로 분리능 Rs ≥ 1.5 확보, "
        "③ 타깃이 없는 매트릭스 블랭크 확보가 어려우므로 최소한 표준액과 시료의 "
        "피크 스펙트럼 일치 여부는 확인할 것."
    )
    parts.append("")

    if simulated:
        parts.append("## 8. (모의 전용) 입력한 참값 대비 회수율")
        parts.append("")
        parts.append(
            "시뮬레이션에 넣어 둔 정답을 분석 코드가 얼마나 되찾는지 확인하는 표다. "
            "첨가법 회수율이 100%에 가깝고 외부 검량선 회수율이 벗어난다면, "
            "코드와 방법론이 의도대로 작동한다는 뜻이다."
        )
        parts.append("")
        parts.append(recovery_vs_truth_table(results))
        parts.append("")

    warnings: list[str] = list(extra_warnings or [])
    for r in results:
        warnings.extend(r.warnings())
    for cal in cals.values():
        if not cal.passes_linearity:
            if cal.source == "manual":
                warnings.append(
                    f"[검량선/{cal.compound.name_ko}] 검량선을 계수로 직접 입력했습니다. "
                    "직선성·LOD/LOQ가 평가되지 않았으므로, 가능하면 CSV에 "
                    "group=calib 행을 넣어 다시 돌리세요."
                )
            else:
                warnings.append(
                    f"[검량선/{cal.compound.name_ko}] R²={cal.r2:.5f} 로 "
                    f"기준({QC.min_r2}) 미달입니다."
                )
    for d in DRINKS.values():
        if not d.label_verified:
            warnings.append(
                f"[{d.name_ko}] 라벨 카페인 표시량이 미확인 상태입니다. "
                "config.py 의 label_caffeine_mg 를 실물 캔 값으로 고치고 "
                "label_verified=True 로 바꾸세요."
            )

    parts.append("## 점검이 필요한 항목")
    parts.append("")
    if warnings:
        for w in dict.fromkeys(warnings):
            parts.append(f"- {w}")
    else:
        parts.append("- 자동 점검에서 걸린 항목이 없습니다.")
    parts.append("")

    if figures:
        parts.append("## 생성된 그림")
        parts.append("")
        for p in figures:
            parts.append(f"- `{p.name}`")
        parts.append("")

    return "\n".join(parts)


def results_dataframe(results: list[StandardAdditionResult]) -> pd.DataFrame:
    """엑셀로 옮기기 좋은 평면 표."""
    rows = []
    for r in results:
        lo, hi = r.ci_c_vial_ppm()
        rows.append({
            "drink_key": r.drink.key,
            "drink": r.drink.name_ko,
            "compound_key": r.compound.key,
            "compound": r.compound.name_ko,
            "sa_slope": r.fit.slope,
            "sa_intercept_A0": r.fit.intercept,
            "sa_r2": r.fit.r2,
            "x_intercept_ppm": r.fit.x_intercept,
            "c_vial_ppm": r.c_vial_ppm,
            "se_c_vial_ppm": r.se_c_vial_ppm,
            "ci_low_ppm": lo,
            "ci_high_ppm": hi,
            "c_drink_ppm": r.c_drink_ppm,
            "mg_per_serving": r.mg_per_serving,
            "se_mg_per_serving": r.se_mg_per_serving,
            "c_apparent_vial_ppm": r.c_apparent_vial_ppm,
            "c_apparent_drink_ppm": r.c_apparent_drink_ppm,
            "conversion_constant": r.conversion_constant,
            "se_conversion_constant": r.se_conversion_constant,
            "external_bias_percent": r.bias_percent,
            "matrix_effect_percent": r.matrix_effect_percent,
            "x_intercept_rsd_percent": r.x_intercept_rsd_percent,
        })
    return pd.DataFrame(rows)


def calibration_dataframe(cals: dict[str, CalibrationResult]) -> pd.DataFrame:
    rows = []
    for key in COMPOUND_ORDER:
        cal = cals.get(key)
        if cal is None:
            continue
        rows.append({
            "compound_key": key,
            "compound": cal.compound.name_ko,
            "source": cal.source,
            "slope": cal.slope,
            "intercept": cal.intercept,
            "r2": cal.r2,
            "se_slope": cal.se_slope,
            "se_intercept": cal.fit.se_intercept if cal.fit else float("nan"),
            "lod_ppm": cal.lod_ppm,
            "loq_ppm": cal.loq_ppm,
            "n_levels": cal.fit.n if cal.fit else 0,
        })
    return pd.DataFrame(rows)
