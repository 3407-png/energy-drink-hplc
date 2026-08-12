"""제로 에너지 드링크 HPLC 동시분석 - 표준물 첨가법 정량 도구.

  hplc.config    실험 조건, 화합물 물성, 시료 정보
  hplc.stats     선형 회귀 및 x절편 오차 전파
  hplc.dataio    피크 면적 CSV 입출력/검증
  hplc.analysis  검량선, 표준물 첨가법, 변환 상수
  hplc.simulate  모의 데이터 생성 (SIMULATION ONLY)
  hplc.plots     그림
  hplc.report    표 및 마크다운 보고서
  hplc.cli       명령줄 진입점
"""

__version__ = "1.0.0"

__all__ = [
    "analysis",
    "config",
    "dataio",
    "plots",
    "report",
    "simulate",
    "stats",
]
