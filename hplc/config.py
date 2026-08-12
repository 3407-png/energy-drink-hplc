"""실험 조건 · 시료 정보 · 화합물 물성 정의.

여기 있는 값은 전부 '연구 계획서에 적힌 조건'을 코드로 옮긴 것이다.
실제 실험에서 조건을 바꿨다면 반드시 이 파일을 먼저 수정할 것.
"""

from __future__ import annotations

from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# 1. HPLC 구동 조건 (연구 목적 - 첫째)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HPLCConditions:
    elution_mode: str = "Isocratic"
    mobile_phase_a: str = "0.1% 인산 수용액"
    mobile_phase_b: str = "아세토나이트릴"
    ratio_a: int = 80
    ratio_b: int = 20
    flow_rate_mL_min: float = 1.0
    oven_temp_C: float = 35.0
    detection_nm: float = 230.0
    column: str = "C18 (4.6 x 250 mm, 5 um)"
    injection_uL: float = 10.0
    run_time_min: float = 10.0

    def summary(self) -> str:
        return (
            f"{self.elution_mode} / {self.mobile_phase_a}:{self.mobile_phase_b} = "
            f"{self.ratio_a}:{self.ratio_b} / {self.flow_rate_mL_min} mL/min / "
            f"{self.oven_temp_C} C / {self.detection_nm} nm"
        )


HPLC = HPLCConditions()


# ---------------------------------------------------------------------------
# 2. 분석 대상 화합물
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Compound:
    key: str
    name_ko: str
    name_en: str
    mw: float                    # g/mol
    pka: float | None            # 산해리상수 (역상 머무름 거동 해석용)
    lambda_max_nm: float         # 최대 흡수 파장
    expected_rt_min: float       # 위 조건에서 예상되는 머무름 시간
    rt_window_min: float = 0.30  # 피크 동정 허용 폭 (+-)

    # ADI: 1일 섭취 허용량 (mg/kg bw/day). 카페인은 ADI 대신 식약처 권고 최대량.
    adi_mg_per_kg: float | None = None
    adi_note: str = ""


CAFFEINE = Compound(
    key="caffeine",
    name_ko="카페인",
    name_en="Caffeine",
    mw=194.19,
    pka=None,             # 중성 (약염기, pKa ~0.7 로 사실상 비이온성)
    lambda_max_nm=273.0,
    expected_rt_min=3.5,
    adi_mg_per_kg=None,
    adi_note="식약처 1일 최대 섭취권고량: 성인 400 mg, 청소년 체중 1 kg당 2.5 mg",
)

SODIUM_BENZOATE = Compound(
    key="sodium_benzoate",
    name_ko="소듐 벤조에이트",
    name_en="Sodium benzoate",
    mw=144.11,
    pka=4.20,             # 벤조산 기준. 이동상 pH~2.3 에서 비해리 -> 머무름 증가
    lambda_max_nm=225.0,
    expected_rt_min=5.5,
    adi_mg_per_kg=5.0,
    adi_note="JECFA/EFSA ADI 0-5 mg/kg bw/day (벤조산 환산)",
)

ACESULFAME_K = Compound(
    key="acesulfame_k",
    name_ko="아세설팜칼륨",
    name_en="Acesulfame K",
    mw=201.24,
    pka=2.0,              # 강산성 -> 역상에서 머무름 짧음
    lambda_max_nm=227.0,
    expected_rt_min=2.2,
    adi_mg_per_kg=15.0,
    adi_note="JECFA ADI 0-15 mg/kg bw/day",
)

COMPOUNDS: dict[str, Compound] = {
    c.key: c for c in (ACESULFAME_K, CAFFEINE, SODIUM_BENZOATE)
}

# 크로마토그램 출력 순서 = 머무름 시간 순
COMPOUND_ORDER = [c.key for c in sorted(COMPOUNDS.values(), key=lambda c: c.expected_rt_min)]


# ---------------------------------------------------------------------------
# 3. 시료 전처리 조건 (탐구 1)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PrepConditions:
    stock_ppm: float = 1000.0           # 표준원액 100 mg / 100 mL
    final_volume_mL: float = 10.0       # 모든 바이알 최종 부피
    drink_aliquot_mL: float = 1.0       # 바이알당 음료 원액 양
    spike_levels_ppm: tuple[float, ...] = (0.0, 20.0, 40.0, 60.0)
    injections_calib: int = 1
    injections_sample: int = 4
    filter_um: float = 0.45

    @property
    def dilution_factor(self) -> float:
        """음료 원액 -> 바이알 희석배수. 1 mL -> 10 mL 이므로 10."""
        return self.final_volume_mL / self.drink_aliquot_mL

    def spike_volume_mL(self, spike_ppm: float) -> float:
        """목표 첨가 농도를 만들기 위해 넣어야 하는 1000 ppm 표준액 부피."""
        return spike_ppm * self.final_volume_mL / self.stock_ppm


PREP = PrepConditions()


# ---------------------------------------------------------------------------
# 4. 시료 (제로 에너지 드링크 3종)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Drink:
    key: str
    name_ko: str
    serving_mL: float
    # 라벨 표기 카페인량 (mg/1회 제공량). 함량 표기 의무 대상.
    label_caffeine_mg: float | None = None
    # 라벨 값을 실물 캔에서 직접 확인했는가?
    label_verified: bool = False
    note: str = ""


# !! 중요 !!
# 아래 라벨 값은 '자리표시자(placeholder)'다. 실험 당일 캔 뒷면을 직접 읽고
# label_caffeine_mg 를 고친 뒤 label_verified=True 로 바꿀 것.
# 검증 안 된 라벨값으로 계산하면 '연구 목적 셋째'(라벨 대조 검증)가 통째로 무의미해진다.
DRINKS: dict[str, Drink] = {
    "monster": Drink(
        key="monster",
        name_ko="몬스터 에너지 울트라",
        serving_mL=355.0,
        label_caffeine_mg=100.0,
        label_verified=False,
        note="라벨 확인 필요",
    ),
    "hotsix": Drink(
        key="hotsix",
        name_ko="핫식스 더킹 크러쉬 피치",
        serving_mL=355.0,
        label_caffeine_mg=100.0,
        label_verified=False,
        note="라벨 확인 필요",
    ),
    "wisely": Drink(
        key="wisely",
        name_ko="와이즐리 에너지 제로 슈가",
        serving_mL=355.0,
        label_caffeine_mg=120.0,
        label_verified=False,
        note="라벨 확인 필요",
    ),
}

# CSV의 group 열에서 검량선 시료를 가리키는 값
CALIB_GROUP = "calib"
SAMPLE_GROUP = "sample"


# ---------------------------------------------------------------------------
# 5. 국내 기준값 (결론/제언 작성용 참고치)
# ---------------------------------------------------------------------------

# 고카페인 함유 표시 대상 기준: 총 카페인 함량이 1 mL당 0.15 mg (=150 ppm) 이상
HIGH_CAFFEINE_THRESHOLD_PPM = 150.0

# 청소년(체중 50 kg 가정) 1일 카페인 권고 상한
TEEN_BODY_WEIGHT_KG = 50.0
CAFFEINE_MG_PER_KG_LIMIT = 2.5


@dataclass(frozen=True)
class Thresholds:
    high_caffeine_ppm: float = HIGH_CAFFEINE_THRESHOLD_PPM
    teen_bw_kg: float = TEEN_BODY_WEIGHT_KG
    caffeine_limit_mg_per_kg: float = CAFFEINE_MG_PER_KG_LIMIT

    @property
    def teen_daily_caffeine_mg(self) -> float:
        return self.teen_bw_kg * self.caffeine_limit_mg_per_kg


THRESHOLDS = Thresholds()


# ---------------------------------------------------------------------------
# 6. 분석 품질 판정 기준 (ICH Q2(R2) 관례값)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class QualityCriteria:
    min_r2: float = 0.995            # 검량선 직선성
    max_rsd_percent: float = 5.0     # 반복 주입 정밀도 (면적 RSD)
    min_resolution: float = 1.5      # 인접 피크 분리능
    # 표준물 첨가법 외삽 신뢰 구간: x절편 상대표준오차가 이 값을 넘으면 경고
    max_x_intercept_rsd_percent: float = 10.0
    confidence_level: float = 0.95


QC = QualityCriteria()
