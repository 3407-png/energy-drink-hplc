"""선형 회귀 및 표준물 첨가법 x절편 오차 전파.

표준물 첨가법의 핵심은 '회귀직선을 데이터 바깥으로 외삽한 지점(x절편)'을 쓴다는
것이다. 외삽은 내삽보다 오차가 크므로, 기울기·절편만 구하고 끝내면 안 되고
x절편의 불확도까지 반드시 같이 보고해야 한다. (ICH Q2(R2) 정확성/정밀성 항목)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import stats as sps


@dataclass(frozen=True)
class LinearFit:
    """최소제곱 직선 y = slope * x + intercept 와 그 통계량."""

    x: np.ndarray
    y: np.ndarray
    slope: float
    intercept: float
    r2: float
    s_yx: float          # 회귀 잔차 표준편차 (standard error of regression)
    se_slope: float
    se_intercept: float
    n: int

    # -- 기본량 ------------------------------------------------------------
    @property
    def dof(self) -> int:
        return self.n - 2

    @property
    def sxx(self) -> float:
        return float(np.sum((self.x - self.x.mean()) ** 2))

    def predict(self, x) -> np.ndarray:
        return self.slope * np.asarray(x, dtype=float) + self.intercept

    @property
    def residuals(self) -> np.ndarray:
        return self.y - self.predict(self.x)

    # -- 표준물 첨가법 x절편 -----------------------------------------------
    @property
    def x_intercept(self) -> float:
        """y=0 이 되는 x. 표준물 첨가법에서 |x절편| = 시료 중 원래 농도."""
        if self.slope == 0:
            return float("nan")
        return -self.intercept / self.slope

    @property
    def se_x_intercept(self) -> float:
        """x절편의 표준오차.

        s_x0 = (s_yx / |m|) * sqrt( 1/n + ybar^2 / (m^2 * Sxx) )

        Miller & Miller, *Statistics and Chemometrics for Analytical
        Chemistry* 의 표준물 첨가법 외삽 오차식.
        """
        if self.slope == 0 or self.dof <= 0:
            return float("nan")
        ybar = float(self.y.mean())
        inside = 1.0 / self.n + (ybar ** 2) / (self.slope ** 2 * self.sxx)
        return (self.s_yx / abs(self.slope)) * math.sqrt(inside)

    def x_intercept_ci(self, level: float = 0.95) -> tuple[float, float]:
        """x절편의 신뢰구간 (기본 95%)."""
        if self.dof <= 0 or not math.isfinite(self.se_x_intercept):
            return (float("nan"), float("nan"))
        t = float(sps.t.ppf(0.5 + level / 2.0, self.dof))
        half = t * self.se_x_intercept
        return (self.x_intercept - half, self.x_intercept + half)

    # -- 외부 검량선용 역산 -------------------------------------------------
    def concentration_from_area(self, area: float) -> float:
        """검량선을 이용해 면적 -> 농도. (A - b) / m"""
        if self.slope == 0:
            return float("nan")
        return (area - self.intercept) / self.slope

    # -- ICH Q2(R2) 검출/정량 한계 -----------------------------------------
    @property
    def lod(self) -> float:
        """검출한계 3.3 * sigma / S. sigma 는 y절편의 표준오차를 사용."""
        if self.slope == 0:
            return float("nan")
        return 3.3 * self.se_intercept / abs(self.slope)

    @property
    def loq(self) -> float:
        """정량한계 10 * sigma / S."""
        if self.slope == 0:
            return float("nan")
        return 10.0 * self.se_intercept / abs(self.slope)


def linear_fit(x, y) -> LinearFit:
    """가중치 없는 최소제곱 직선 적합."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError(f"x, y 길이가 다릅니다: {x.shape} vs {y.shape}")
    n = x.size
    if n < 3:
        raise ValueError(
            f"회귀에 최소 3점이 필요합니다 (현재 {n}점). "
            "표준물 첨가 농도 수준을 확인하세요."
        )
    if np.allclose(x, x[0]):
        raise ValueError("모든 x값이 동일합니다 - 첨가 농도가 변하지 않았습니다.")

    xbar, ybar = x.mean(), y.mean()
    sxx = float(np.sum((x - xbar) ** 2))
    sxy = float(np.sum((x - xbar) * (y - ybar)))
    syy = float(np.sum((y - ybar) ** 2))

    slope = sxy / sxx
    intercept = ybar - slope * xbar

    ss_res = float(np.sum((y - (slope * x + intercept)) ** 2))
    r2 = 1.0 - ss_res / syy if syy > 0 else float("nan")

    dof = n - 2
    s_yx = math.sqrt(ss_res / dof) if dof > 0 else float("nan")
    se_slope = s_yx / math.sqrt(sxx) if dof > 0 else float("nan")
    se_intercept = (
        s_yx * math.sqrt(1.0 / n + xbar ** 2 / sxx) if dof > 0 else float("nan")
    )

    return LinearFit(
        x=x, y=y,
        slope=slope, intercept=intercept, r2=r2,
        s_yx=s_yx, se_slope=se_slope, se_intercept=se_intercept, n=n,
    )


def rsd_percent(values) -> float:
    """상대표준편차(%) = 표준편차 / 평균 * 100. 반복 주입 정밀도 평가용."""
    v = np.asarray(values, dtype=float)
    if v.size < 2:
        return float("nan")
    m = v.mean()
    if m == 0:
        return float("nan")
    return float(v.std(ddof=1) / m * 100.0)


def resolution(rt1: float, w1: float, rt2: float, w2: float) -> float:
    """인접 두 피크의 분리능 Rs = 2(t2 - t1) / (w1 + w2). w 는 베이스라인 폭."""
    if (w1 + w2) == 0:
        return float("nan")
    return 2.0 * abs(rt2 - rt1) / (w1 + w2)


def propagate_ratio(
    a: float, sa: float, b: float, sb: float
) -> tuple[float, float]:
    """비율 a/b 와 그 표준오차 (독립 가정, 상대오차 제곱합)."""
    if b == 0:
        return (float("nan"), float("nan"))
    ratio = a / b
    if a == 0:
        return (0.0, float("nan"))
    rel = math.sqrt((sa / a) ** 2 + (sb / b) ** 2)
    return (ratio, abs(ratio) * rel)
