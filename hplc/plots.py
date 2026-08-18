"""그림 생성 (보고서 [그림 7] 이후에 넣을 것들).

한글 폰트가 설치되어 있으면 한글 라벨을, 없으면 자동으로 영문 라벨을 쓴다.
(matplotlib 은 폰트가 없으면 글자를 네모(tofu)로 그려 버리기 때문에,
 보고서에 그대로 붙였다가 낭패 보는 걸 막으려는 처리다.)
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

from .analysis import CalibrationResult, StandardAdditionResult
from .config import COMPOUND_ORDER, COMPOUNDS, HPLC

# ---------------------------------------------------------------------------
# 폰트
# ---------------------------------------------------------------------------

_KOREAN_FONT_CANDIDATES = [
    "NanumGothic", "NanumBarunGothic", "Malgun Gothic", "AppleGothic",
    "Noto Sans CJK KR", "Noto Sans KR", "Source Han Sans KR", "UnDotum",
]

_HAS_KOREAN = False


def setup_style() -> bool:
    """한글 폰트를 찾아 설정하고, 성공 여부를 돌려준다."""
    global _HAS_KOREAN
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in _KOREAN_FONT_CANDIDATES:
        if name in available:
            plt.rcParams["font.family"] = name
            _HAS_KOREAN = True
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams.update(
        {
            "figure.dpi": 130,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": "--",
            "axes.axisbelow": True,
            "font.size": 9,
            "axes.titlesize": 10,
            "legend.frameon": False,
        }
    )
    return _HAS_KOREAN


def L(ko: str, en: str) -> str:
    """한글 폰트가 있으면 ko, 없으면 en."""
    return ko if _HAS_KOREAN else en


def _cname(cmp_key: str) -> str:
    c = COMPOUNDS[cmp_key]
    return L(c.name_ko, c.name_en)


# 성분별 고정 색 (모든 그림에서 동일하게 유지)
# acesulfame_k 는 현재 분석 대상이 아니지만, 다시 넣을 때 색이 바뀌지 않도록 남겨 둔다.
COLORS = {
    "acesulfame_k": "#2E7D8F",
    "caffeine": "#C05746",
    "sodium_benzoate": "#6B7F3E",
}
DRINK_COLORS = ["#2E7D8F", "#C05746", "#6B7F3E"]


# ---------------------------------------------------------------------------
# 1. 외부 검량선
# ---------------------------------------------------------------------------

def plot_calibrations(
    calibrations: dict[str, CalibrationResult], outdir: Path
) -> Path:
    keys = [k for k in COMPOUND_ORDER if k in calibrations]
    fig, axes = plt.subplots(1, len(keys), figsize=(4.0 * len(keys), 3.4))
    if len(keys) == 1:
        axes = [axes]

    for ax, key in zip(axes, keys):
        cal = calibrations[key]
        fit = cal.fit
        color = COLORS[key]
        ax.scatter(fit.x, fit.y, s=34, color=color, zorder=3,
                   label=L("실측", "measured"))
        xs = np.linspace(0, max(fit.x) * 1.05, 100)
        ax.plot(xs, fit.predict(xs), color=color, lw=1.4, alpha=0.85)
        ax.set_title(_cname(key))
        ax.set_xlabel(L("농도 (ppm)", "Concentration (ppm)"))
        ax.set_ylabel(L("피크 면적", "Peak area"))
        ax.text(
            0.04, 0.95,
            f"y = {fit.slope:,.0f}x + {fit.intercept:,.0f}\n"
            f"$R^2$ = {fit.r2:.5f}\n"
            f"LOD = {cal.lod_ppm:.2f} ppm\nLOQ = {cal.loq_ppm:.2f} ppm",
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
        )
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    fig.suptitle(
        L(f"외부 검량선 (증류수 바탕, {HPLC.detection_nm:.0f} nm)",
          f"External calibration (aqueous, {HPLC.detection_nm:.0f} nm)"),
        y=1.02,
    )
    path = outdir / "fig_calibration_curves.png"
    fig.savefig(path)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 2. 표준물 첨가법 (x절편 외삽)
# ---------------------------------------------------------------------------

def plot_standard_addition(
    results: list[StandardAdditionResult], outdir: Path
) -> list[Path]:
    paths: list[Path] = []
    by_drink: dict[str, list[StandardAdditionResult]] = {}
    for r in results:
        by_drink.setdefault(r.drink.key, []).append(r)

    for drink_key, rs in by_drink.items():
        rs = sorted(rs, key=lambda r: r.compound.expected_rt_min)
        fig, axes = plt.subplots(1, len(rs), figsize=(4.0 * len(rs), 3.6))
        if len(rs) == 1:
            axes = [axes]

        for ax, r in zip(axes, rs):
            fit = r.fit
            color = COLORS[r.compound.key]
            xint = fit.x_intercept

            ax.scatter(fit.x, fit.y, s=36, color=color, zorder=4,
                       label=L("실측 (첨가 시료)", "measured"))
            # 실측 구간은 실선, 데이터 바깥으로 나가는 외삽 구간은 점선
            xs = np.linspace(0, max(fit.x) * 1.05, 100)
            ax.plot(xs, fit.predict(xs), color=color, lw=1.4, zorder=3,
                    label=L("회귀직선", "regression"))
            xs_ext = np.linspace(min(xint * 1.25, -1.0), 0, 80)
            ax.plot(xs_ext, fit.predict(xs_ext), color=color, lw=1.4,
                    ls="--", dashes=(4, 3), alpha=0.9, zorder=3,
                    label=L("외삽", "extrapolation"))

            ax.axhline(0, color="0.35", lw=0.8)
            ax.axvline(0, color="0.35", lw=0.8)
            ax.scatter([xint], [0], marker="v", s=70, color="black", zorder=5)

            lo, hi = r.ci_c_vial_ppm()
            if math.isfinite(lo) and math.isfinite(hi):
                ax.axvspan(-hi, -lo, color=color, alpha=0.12, lw=0)

            ax.annotate(
                L(f"x절편 = {xint:.2f} ppm\n→ 시료 농도 {r.c_vial_ppm:.2f} ppm",
                  f"x-int = {xint:.2f} ppm\n→ {r.c_vial_ppm:.2f} ppm"),
                xy=(xint, 0), xytext=(0.05, 0.72), textcoords="axes fraction",
                fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.8, color="black"),
            )
            ax.text(
                0.05, 0.95,
                f"y = {fit.slope:,.0f}x + {fit.intercept:,.0f}\n$R^2$ = {fit.r2:.5f}",
                transform=ax.transAxes, va="top", fontsize=8,
            )
            ax.set_title(_cname(r.compound.key))
            ax.set_xlabel(L("첨가 농도 (ppm)", "Added concentration (ppm)"))
            ax.set_ylabel(L("피크 면적", "Peak area"))
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

        title_drink = L(rs[0].drink.name_ko, rs[0].drink.key)
        fig.suptitle(
            L(f"표준물 첨가법 - {title_drink}",
              f"Standard addition - {title_drink}"),
            y=1.02,
        )
        path = outdir / f"fig_standard_addition_{drink_key}.png"
        fig.savefig(path)
        plt.close(fig)
        paths.append(path)
    return paths


# ---------------------------------------------------------------------------
# 3. 변환 상수 / 매트릭스 효과
# ---------------------------------------------------------------------------

def plot_conversion_constants(
    results: list[StandardAdditionResult], outdir: Path
) -> Path:
    drinks = sorted({r.drink.key for r in results})
    keys = [k for k in COMPOUND_ORDER if any(r.compound.key == k for r in results)]
    lookup = {(r.drink.key, r.compound.key): r for r in results}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.8))
    width = 0.8 / max(len(drinks), 1)
    idx = np.arange(len(keys))

    from .config import DRINKS

    for i, dk in enumerate(drinks):
        f_vals, e_vals, m_vals = [], [], []
        for k in keys:
            r = lookup.get((dk, k))
            f_vals.append(r.conversion_constant if r else np.nan)
            e_vals.append(r.se_conversion_constant if r else np.nan)
            m_vals.append(r.matrix_effect_percent if r else np.nan)
        label = L(DRINKS[dk].name_ko, dk) if dk in DRINKS else dk
        color = DRINK_COLORS[i % len(DRINK_COLORS)]
        pos = idx + i * width - 0.4 + width / 2
        # 변환 상수는 1 근처의 비율이므로 막대(0 기준)보다 1 기준 편차로 보여준다
        ax1.vlines(pos, 1.0, f_vals, color=color, lw=1.2, alpha=0.8)
        ax1.errorbar(pos, f_vals, yerr=e_vals, fmt="o", ms=6, capsize=3,
                     color=color, ecolor=color, elinewidth=1.0, label=label)
        ax2.bar(pos, m_vals, width * 0.9, color=color, alpha=0.85, label=label)

    ax1.axhline(1.0, color="0.3", lw=1.0, ls="--")
    ax1.margins(y=0.35)
    ax1.set_xticks(idx, [_cname(k) for k in keys])
    ax1.set_ylabel(L("변환 상수 f", "Conversion constant f"))
    ax1.set_title(L("변환 상수 (참농도/겉보기농도)",
                    "Conversion constant (true/apparent)"))
    ax1.legend(fontsize=8)

    ax2.axhline(0.0, color="0.3", lw=1.0, ls="--")
    ax2.set_xticks(idx, [_cname(k) for k in keys])
    ax2.set_ylabel(L("기울기 변화율 (%)", "Slope change (%)"))
    ax2.set_title(L("곱셈형 매트릭스 효과", "Multiplicative matrix effect"))

    path = outdir / "fig_conversion_constants.png"
    fig.savefig(path)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 4. 크로마토그램 (모의)
# ---------------------------------------------------------------------------

def plot_chromatogram(
    t: np.ndarray, signal: np.ndarray, peaks: list[dict], title: str, outdir: Path,
    filename: str, simulated: bool = True,
) -> Path:
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    ax.plot(t, signal, color="#1f3b4d", lw=0.9)

    top = float(signal.max())
    ax.set_ylim(-0.05 * top, top * 1.30)

    for p in peaks:
        i = int(np.argmin(np.abs(t - p["rt"])))
        name = L(p["name"], p.get("name_en", ""))
        if p.get("target"):
            ax.annotate(
                f"{name}\n{p['rt']:.2f} min",
                xy=(t[i], signal[i]), xytext=(0, 10), textcoords="offset points",
                ha="center", va="bottom", fontsize=8,
            )
        else:
            ax.annotate(
                name, xy=(t[i], signal[i]), xytext=(0, 7),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=7, color="0.45",
            )

    ax.set_xlabel(L("머무름 시간 (분)", "Retention time (min)"))
    ax.set_ylabel(L("흡광도 (mAU)", "Absorbance (mAU)"))
    ax.set_title(title, pad=12)
    ax.set_xlim(0, t.max())
    if simulated:
        ax.text(
            0.99, 0.95, "SIMULATED", transform=ax.transAxes, ha="right", va="top",
            fontsize=16, color="red", alpha=0.28, weight="bold",
        )
    path = outdir / filename
    fig.savefig(path)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 5. 잔차 (직선성 진단)
# ---------------------------------------------------------------------------

def plot_residuals(results: list[StandardAdditionResult], outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    markers = ["o", "s", "^"]
    drinks = sorted({r.drink.key for r in results})
    for r in results:
        di = drinks.index(r.drink.key)
        rel = r.fit.residuals / np.maximum(np.abs(r.fit.y), 1) * 100
        ax.scatter(
            r.fit.x, rel, marker=markers[di % len(markers)],
            color=COLORS[r.compound.key], s=30, alpha=0.85,
        )
    ax.axhline(0, color="0.3", lw=1.0)
    ax.set_xlabel(L("첨가 농도 (ppm)", "Added concentration (ppm)"))
    ax.set_ylabel(L("상대 잔차 (%)", "Relative residual (%)"))
    ax.set_title(L("표준물 첨가법 회귀 잔차 (직선성 진단)",
                   "Standard addition residuals"))
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=COLORS[k], label=_cname(k))
        for k in COMPOUND_ORDER if any(r.compound.key == k for r in results)
    ]
    ax.legend(handles=handles, fontsize=8)
    path = outdir / "fig_residuals.png"
    fig.savefig(path)
    plt.close(fig)
    return path
