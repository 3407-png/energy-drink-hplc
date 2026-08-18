"""모의 데이터 생성기 (SIMULATION ONLY).

===========================================================================
 이 모듈이 만드는 숫자는 실측값이 아니다.
 목적은 두 가지뿐이다.
   (1) 실험 전에 분석 코드가 제대로 도는지 검증한다.
   (2) "이 조건이면 대략 이런 크로마토그램과 결과가 나온다"는 예상치를 잡는다.
 여기서 나온 값을 연구 결과(7장)에 실측값으로 적으면 데이터 조작이다.
 보고서에 쓸 거라면 반드시 '예상 결과(모의 계산)'로 명시할 것.
===========================================================================

모형
    검량선 바이알 :  A = m * C_add + b_cal + noise
    음료   바이알 :  A = m * f_resp * (C_native + C_add) + b_off + noise

    m       순수 용매에서의 응답계수 (면적/ppm)
    f_resp  곱셈형 매트릭스 효과. 매트릭스가 단위 농도당 응답을 바꾸는 정도.
            -> 표준물 첨가법이 보정해 주는 성분.
    b_off   덧셈형 간섭. 같은 머무름 시간대에 겹쳐 나오는 성분(리보플라빈 등)이
            더해 놓은 면적.
            -> 표준물 첨가법으로 보정되지 '않는' 성분. x절편을 부풀린다.

아래 상수들은 문헌에 보고된 일반적인 범위를 참고해 임의로 잡은 가정값이다.
실제 제품의 함량이 아니다.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .config import (
    CALIB_GROUP,
    COMPOUND_ORDER,
    COMPOUNDS,
    DRINKS,
    HPLC,
    PREP,
    SAMPLE_GROUP,
    VOID_TIME_MIN,
)

SIM_BANNER = [
    "# ===================================================================",
    "# SIMULATED DATA - NOT EXPERIMENTAL MEASUREMENT",
    "# 모의 생성 데이터입니다. 실측값이 아닙니다.",
    "# 코드 검증 및 예상 결과 확인용으로만 사용하십시오.",
    "# ===================================================================",
]


# ---------------------------------------------------------------------------
# 가정값
# ---------------------------------------------------------------------------

# 순수 용매에서의 응답계수 (면적 단위 / ppm). 230 nm 흡광 세기 차이를 반영.
RESPONSE_FACTOR: dict[str, float] = {
    "caffeine": 21_000.0,
    "sodium_benzoate": 34_000.0,
}

# 검량선 y절편 (적분 바탕선 오프셋). 이상적으로는 0 근처.
CALIB_INTERCEPT: dict[str, float] = {
    "caffeine": 900.0,
    "sodium_benzoate": 1_500.0,
}

# 음료 원액 중 '참' 농도 (ppm = mg/L) — 가정값
TRUE_DRINK_PPM: dict[str, dict[str, float]] = {
    "monster":  {"caffeine": 282.0, "sodium_benzoate": 250.0},
    "netflix":  {"caffeine": 300.0, "sodium_benzoate": 220.0},
    "wisely":   {"caffeine": 340.0, "sodium_benzoate": 180.0},
}

# 곱셈형 매트릭스 효과 (1.0 = 영향 없음)
MATRIX_RESPONSE: dict[str, dict[str, float]] = {
    "monster":  {"caffeine": 0.93, "sodium_benzoate": 0.98},
    "netflix":  {"caffeine": 0.89, "sodium_benzoate": 0.95},
    "wisely":   {"caffeine": 0.97, "sodium_benzoate": 1.02},
}

# 덧셈형 공용리 간섭이 더해 놓는 면적 (리보플라빈/유기산 등)
ADDITIVE_OFFSET: dict[str, dict[str, float]] = {
    "monster":  {"caffeine": 12_000.0, "sodium_benzoate": 3_000.0},
    "netflix":  {"caffeine": 22_000.0, "sodium_benzoate": 4_500.0},
    "wisely":   {"caffeine": 6_000.0,  "sodium_benzoate": 2_000.0},
}

# 잡음
AREA_RSD = 0.012          # 주입 재현성 1.2%
BASELINE_NOISE_AREA = 400.0
RT_JITTER_MIN = 0.02


@dataclass
class SimulationSettings:
    seed: int = 20260812
    area_rsd: float = AREA_RSD
    baseline_noise: float = BASELINE_NOISE_AREA
    include_additive: bool = True
    include_matrix: bool = True


# ---------------------------------------------------------------------------
# 피크 면적 데이터
# ---------------------------------------------------------------------------

def simulate_peak_areas(settings: SimulationSettings | None = None) -> pd.DataFrame:
    """탐구 1~2를 그대로 따라간 가상의 피크 면적표를 만든다."""
    st = settings or SimulationSettings()
    rng = np.random.default_rng(st.seed)
    rows: list[dict] = []

    def noisy(area: float) -> float:
        val = area * (1.0 + rng.normal(0.0, st.area_rsd))
        val += rng.normal(0.0, st.baseline_noise)
        return max(val, 0.0)

    # --- 검량선 바이알 (증류수 바탕) ---
    for cmp_key in COMPOUND_ORDER:
        m = RESPONSE_FACTOR[cmp_key]
        b = CALIB_INTERCEPT[cmp_key]
        rt0 = COMPOUNDS[cmp_key].expected_rt_min
        for spike in PREP.spike_levels_ppm:
            for inj in range(1, PREP.injections_calib + 1):
                area = noisy(m * spike + b) if spike > 0 else max(
                    rng.normal(b, st.baseline_noise), 0.0
                )
                rows.append(
                    {
                        "group": CALIB_GROUP,
                        "sample": "STD",
                        "compound": cmp_key,
                        "spike_ppm": spike,
                        "injection": inj,
                        "retention_min": round(
                            rt0 + rng.normal(0, RT_JITTER_MIN), 3
                        ),
                        "peak_area": round(area, 1),
                    }
                )

    # --- 음료 바이알 (표준물 첨가) ---
    for drink_key in DRINKS:
        for cmp_key in COMPOUND_ORDER:
            m = RESPONSE_FACTOR[cmp_key]
            rt0 = COMPOUNDS[cmp_key].expected_rt_min
            c_native_vial = TRUE_DRINK_PPM[drink_key][cmp_key] / PREP.dilution_factor
            f_resp = MATRIX_RESPONSE[drink_key][cmp_key] if st.include_matrix else 1.0
            b_off = ADDITIVE_OFFSET[drink_key][cmp_key] if st.include_additive else 0.0

            for spike in PREP.spike_levels_ppm:
                true_area = m * f_resp * (c_native_vial + spike) + b_off
                for inj in range(1, PREP.injections_sample + 1):
                    rows.append(
                        {
                            "group": SAMPLE_GROUP,
                            "sample": drink_key,
                            "compound": cmp_key,
                            "spike_ppm": spike,
                            "injection": inj,
                            "retention_min": round(
                                rt0 + rng.normal(0, RT_JITTER_MIN), 3
                            ),
                            "peak_area": round(noisy(true_area), 1),
                        }
                    )

    return pd.DataFrame(rows)


def write_simulated_csv(path: str | Path, settings: SimulationSettings | None = None) -> Path:
    """모의 데이터를 SIMULATED 배너가 박힌 CSV로 저장한다."""
    df = simulate_peak_areas(settings)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    st = settings or SimulationSettings()
    with path.open("w", encoding="utf-8-sig") as fh:
        for line in SIM_BANNER:
            fh.write(line + "\n")
        fh.write(f"# seed={st.seed}, matrix={st.include_matrix}, additive={st.include_additive}\n")
        df.to_csv(fh, index=False, lineterminator="\n")
    return path


def truth_table() -> pd.DataFrame:
    """시뮬레이션에 넣은 '정답'. 역산 결과가 이걸 얼마나 되찾는지 보는 용도."""
    rows = []
    for drink_key, d in DRINKS.items():
        for cmp_key in COMPOUND_ORDER:
            c_drink = TRUE_DRINK_PPM[drink_key][cmp_key]
            rows.append(
                {
                    "drink": drink_key,
                    "drink_name": d.name_ko,
                    "compound": cmp_key,
                    "compound_name": COMPOUNDS[cmp_key].name_ko,
                    "true_drink_ppm": c_drink,
                    "true_vial_ppm": c_drink / PREP.dilution_factor,
                    "true_mg_per_serving": c_drink * d.serving_mL / 1000.0,
                    "matrix_response": MATRIX_RESPONSE[drink_key][cmp_key],
                    "additive_offset_area": ADDITIVE_OFFSET[drink_key][cmp_key],
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 크로마토그램 (그림용)
# ---------------------------------------------------------------------------

# 면적(임의 단위) -> 검출기 응답(mAU) 환산 계수.
# 40 ppm 표준액의 카페인 피크 높이가 대략 300 mAU 가 되도록 맞춘 값이다.
AREA_TO_MAU = 5.0e-5


def _gaussian(t: np.ndarray, center: float, area: float, width: float) -> np.ndarray:
    """면적이 주어진 가우시안 피크를 mAU 단위 높이로 그린다."""
    sigma = width / 2.355  # FWHM -> sigma
    amplitude = area * AREA_TO_MAU / (sigma * np.sqrt(2 * np.pi))
    return amplitude * np.exp(-0.5 * ((t - center) / sigma) ** 2)


def simulate_chromatogram(
    drink_key: str | None = None,
    spike_ppm: float = 0.0,
    settings: SimulationSettings | None = None,
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """가상 크로마토그램 (시간, 신호, 피크 목록)을 만든다.

    drink_key=None 이면 검량선 표준액(증류수 바탕) 크로마토그램.
    """
    st = settings or SimulationSettings()
    rng = np.random.default_rng(st.seed + int(spike_ppm))
    t = np.linspace(0.0, HPLC.run_time_min, 6000)
    signal = np.zeros_like(t)
    peaks: list[dict] = []

    # 용매 피크 (컬럼 공극 시간)
    t0 = VOID_TIME_MIN
    signal += _gaussian(t, t0, 9_000.0, 0.20)
    peaks.append({"name": "용매 피크", "name_en": "solvent front", "rt": t0,
                  "target": False})

    for cmp_key in COMPOUND_ORDER:
        c = COMPOUNDS[cmp_key]
        m = RESPONSE_FACTOR[cmp_key]
        if drink_key is None:
            area = m * spike_ppm + CALIB_INTERCEPT[cmp_key]
        else:
            c_native_vial = TRUE_DRINK_PPM[drink_key][cmp_key] / PREP.dilution_factor
            f_resp = MATRIX_RESPONSE[drink_key][cmp_key] if st.include_matrix else 1.0
            area = m * f_resp * (c_native_vial + spike_ppm)
        width = 0.10 + 0.012 * c.expected_rt_min
        signal += _gaussian(t, c.expected_rt_min, area, width)
        peaks.append({
            "name": c.name_ko, "name_en": c.name_en,
            "rt": c.expected_rt_min, "target": True,
        })

    # 매트릭스 유래 방해 피크 (음료 시료에만)
    if drink_key is not None and st.include_additive:
        # 70:30 조건에서는 매트릭스 성분도 전부 앞으로 당겨진다.
        # 리보플라빈은 카페인(3.8분) 바로 앞에 붙어 공용리 위험이 가장 큰 성분이다.
        matrix_peaks = [
            ("유기산류", "organic acids", 2.7, 30_000.0, 0.22),
            ("리보플라빈", "riboflavin", 3.5, 18_000.0, 0.16),
            ("향료 성분", "flavor cmpd.", 6.8, 12_000.0, 0.28),
        ]
        scale = {"monster": 1.0, "netflix": 1.7, "wisely": 0.45}.get(drink_key, 1.0)
        for name, name_en, rt, area, width in matrix_peaks:
            signal += _gaussian(t, rt, area * scale, width)
            peaks.append({"name": name, "name_en": name_en, "rt": rt,
                          "target": False})

    # 바탕선 드리프트 + 잡음 (mAU)
    signal += 1.0 + 0.6 * t
    signal += rng.normal(0.0, 0.15, size=t.size)
    return t, signal, peaks
