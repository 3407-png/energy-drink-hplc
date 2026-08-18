"""피크 면적 데이터 입출력 및 검증.

HPLC 소프트웨어에서 내보낸 결과를 아래 형식의 CSV 한 장으로 정리해서 넣으면 된다.

    group,sample,compound,spike_ppm,injection,retention_min,peak_area

  group        : calib(검량선용, 증류수 바탕) 또는 sample(음료 바탕)
  sample       : calib 이면 STD, sample 이면 monster / netflix / wisely
  compound     : caffeine / sodium_benzoate
  spike_ppm    : 그 바이알에 첨가한 표준물질 농도 (0, 20, 40, 60)
  injection    : 같은 바이알의 몇 번째 주입인지 (1부터)
  retention_min: 실측 머무름 시간 (분). 모르면 비워둬도 됨
  peak_area    : 피크 면적 (기기 출력 단위 그대로)

피크가 검출되지 않은 경우 peak_area 에 0 을 쓴다(빈칸 아님).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import CALIB_GROUP, COMPOUNDS, DRINKS, PREP, SAMPLE_GROUP

COLUMNS = [
    "group",
    "sample",
    "compound",
    "spike_ppm",
    "injection",
    "retention_min",
    "peak_area",
]


class DataValidationError(ValueError):
    """입력 CSV가 분석 가능한 형태가 아닐 때."""


def write_template(path: str | Path) -> Path:
    """실험 당일 그대로 채워 넣을 빈 CSV 템플릿을 만든다."""
    rows = []
    for cmp_key in COMPOUNDS:
        for spike in PREP.spike_levels_ppm:
            for inj in range(1, PREP.injections_calib + 1):
                rows.append(
                    {
                        "group": CALIB_GROUP,
                        "sample": "STD",
                        "compound": cmp_key,
                        "spike_ppm": spike,
                        "injection": inj,
                        "retention_min": "",
                        "peak_area": "",
                    }
                )
    for drink_key in DRINKS:
        for cmp_key in COMPOUNDS:
            for spike in PREP.spike_levels_ppm:
                for inj in range(1, PREP.injections_sample + 1):
                    rows.append(
                        {
                            "group": SAMPLE_GROUP,
                            "sample": drink_key,
                            "compound": cmp_key,
                            "spike_ppm": spike,
                            "injection": inj,
                            "retention_min": "",
                            "peak_area": "",
                        }
                    )

    df = pd.DataFrame(rows, columns=COLUMNS)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def is_simulated(path: str | Path) -> bool:
    """파일 머리말에 SIMULATED 배너가 있으면 모의 데이터로 본다."""
    path = Path(path)
    if not path.exists():
        return False
    with path.open("r", encoding="utf-8-sig") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "SIMULATED" in line.upper():
                return True
    return False


def load_peak_areas(path: str | Path) -> pd.DataFrame:
    """CSV를 읽고 분석 전에 형식을 검증한다.

    '#' 로 시작하는 머리말 줄은 주석으로 무시한다(모의 데이터 배너용).
    """
    path = Path(path)
    if not path.exists():
        raise DataValidationError(f"파일이 없습니다: {path}")

    df = pd.read_csv(path, comment="#")
    df.columns = [c.strip().lower() for c in df.columns]

    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise DataValidationError(
            f"필수 열이 없습니다: {missing}\n필요한 열: {COLUMNS}"
        )

    for col in ("group", "sample", "compound"):
        df[col] = df[col].astype(str).str.strip()
    for col in ("spike_ppm", "peak_area", "retention_min"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["injection"] = pd.to_numeric(df["injection"], errors="coerce").fillna(1).astype(int)

    # 아직 안 채운 행은 조용히 버리되, 몇 줄인지는 알려준다.
    blank = df["peak_area"].isna()
    if blank.any():
        df = df.loc[~blank].copy()

    if df.empty:
        raise DataValidationError(
            "peak_area 가 채워진 행이 하나도 없습니다. 템플릿을 먼저 채우세요."
        )

    bad_group = set(df["group"]) - {CALIB_GROUP, SAMPLE_GROUP}
    if bad_group:
        raise DataValidationError(
            f"group 열에 알 수 없는 값: {bad_group} (허용: {CALIB_GROUP}, {SAMPLE_GROUP})"
        )

    bad_cmp = set(df["compound"]) - set(COMPOUNDS)
    if bad_cmp:
        raise DataValidationError(
            f"compound 열에 알 수 없는 값: {bad_cmp} (허용: {list(COMPOUNDS)})"
        )

    samples = set(df.loc[df["group"] == SAMPLE_GROUP, "sample"])
    bad_sample = samples - set(DRINKS)
    if bad_sample:
        raise DataValidationError(
            f"sample 열에 알 수 없는 음료: {bad_sample} (허용: {list(DRINKS)})\n"
            "config.py 의 DRINKS 에 추가하거나 CSV의 이름을 맞추세요."
        )

    if (df["peak_area"] < 0).any():
        raise DataValidationError("peak_area 에 음수가 있습니다. 적분 결과를 확인하세요.")

    return df.reset_index(drop=True)


def mean_areas(df: pd.DataFrame) -> pd.DataFrame:
    """반복 주입을 평균으로 묶고 정밀도(RSD%)를 함께 계산한다."""
    from .stats import rsd_percent

    grouped = df.groupby(["group", "sample", "compound", "spike_ppm"], as_index=False)
    out = grouped.agg(
        peak_area_mean=("peak_area", "mean"),
        peak_area_sd=("peak_area", lambda s: s.std(ddof=1)),
        n_injections=("peak_area", "size"),
        retention_min=("retention_min", "mean"),
    )
    rsd = grouped["peak_area"].apply(lambda s: rsd_percent(s.values))
    out["rsd_percent"] = rsd.iloc[:, -1].values if hasattr(rsd, "iloc") else rsd
    return out


def check_retention_times(df: pd.DataFrame) -> list[str]:
    """실측 머무름 시간이 예상 창을 벗어나면 경고 문자열로 돌려준다."""
    warnings: list[str] = []
    if df["retention_min"].isna().all():
        return ["머무름 시간이 입력되지 않아 피크 동정 검증을 건너뜁니다."]

    for cmp_key, sub in df.groupby("compound"):
        c = COMPOUNDS[cmp_key]
        rts = sub["retention_min"].dropna()
        if rts.empty:
            continue
        lo, hi = rts.min(), rts.max()
        if hi - lo > c.rt_window_min * 2:
            warnings.append(
                f"[{c.name_ko}] 머무름 시간 변동이 큽니다 "
                f"({lo:.2f}~{hi:.2f} 분). 컬럼 평형 또는 피크 오동정 확인 필요."
            )
    return warnings
