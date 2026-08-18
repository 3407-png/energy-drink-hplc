"""명령줄 진입점.

    python -m hplc template                 빈 입력 CSV 템플릿 만들기
    python -m hplc design-check             실험 전 첨가 농도 설계 점검
    python -m hplc simulate                 모의 데이터 CSV 만들기 (검증용)
    python -m hplc chromatogram             모의 크로마토그램 그림 만들기
    python -m hplc analyze <csv>            실측/모의 CSV 분석 -> 보고서 + 그림
    python -m hplc demo                     모의 데이터 생성부터 분석까지 한 번에
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import analysis as an
from . import dataio, plots, report
from .config import COMPOUND_ORDER, COMPOUNDS, DRINKS, HPLC, PREP, VOID_TIME_MIN

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA = ROOT / "data"
DEFAULT_OUT = ROOT / "output"


# ---------------------------------------------------------------------------

def cmd_template(args) -> int:
    path = dataio.write_template(args.output)
    print(f"입력 템플릿을 만들었습니다: {path}")
    print()
    print("작성 방법")
    print("  1. HPLC 소프트웨어에서 피크 면적을 내보냅니다.")
    print("  2. peak_area 열을 채웁니다. 피크가 안 나온 경우 0 을 적습니다.")
    print("  3. retention_min 열도 채우면 피크 오동정을 자동으로 점검합니다.")
    print(f"  4. python -m hplc analyze {path}")
    return 0


def cmd_design_check(args) -> int:
    print("=" * 68)
    print("실험 설계 점검 (실험 전에 확인)")
    print("=" * 68)
    print(f"HPLC 조건 : {HPLC.summary()}")
    print(f"컬럼      : {HPLC.column}")
    print()

    print("[1] 예상 머무름 시간 및 분리")
    print(f"    컬럼 공극 시간 t0 = {VOID_TIME_MIN:.1f} 분 (추정)")
    ordered = [COMPOUNDS[k] for k in COMPOUND_ORDER]
    for c in ordered:
        pka = f"pKa {c.pka}" if c.pka is not None else "비이온성"
        k_prime = (c.expected_rt_min - VOID_TIME_MIN) / VOID_TIME_MIN
        print(
            f"  {c.name_ko:<12} RT ~{c.expected_rt_min:>4.1f} 분  "
            f"k' {k_prime:>4.1f}  λmax {c.lambda_max_nm:>5.0f} nm  ({pka})"
        )
        if k_prime < 1.0:
            print(
                "               ⚠ k' < 1 — 컬럼에 거의 붙잡히지 않습니다. "
                "공극 근처에서 나오는 매트릭스 성분과 겹칠 위험이 큽니다."
            )
    for a, b in zip(ordered, ordered[1:]):
        gap = b.expected_rt_min - a.expected_rt_min
        flag = "여유 있음" if gap >= 0.8 else "겹칠 위험 — 분리능 확인 필요"
        print(f"  · {a.name_ko} ↔ {b.name_ko}: 간격 {gap:.1f} 분 ({flag})")
    print()
    print(f"  참고: 검출 파장 {HPLC.detection_nm:.0f} nm 는 두 성분의 절충값입니다.")
    print("        카페인 λmax 는 273 nm 이므로 230 nm 에서는 감도가 손해입니다.")
    print("        분석 성분이 2종으로 줄었으므로, PDA 를 쓸 수 있다면")
    print("        카페인 273 nm / 소듐벤조에이트 225 nm 로 따로 뽑는 편이 낫습니다.")
    print()
    print("  ⚠ 위 머무름 시간은 이동상 70:30 조건의 추정값입니다.")
    print("     첫 크로마토그램을 받으면 config.py 의 expected_rt_min 을")
    print("     실측값으로 바꾸세요. 피크 동정 점검 기준이 됩니다.")
    print()

    print("[2] 표준물 첨가 농도 적정성")
    print(f"    첨가 수준: {', '.join(f'{s:.0f}' for s in PREP.spike_levels_ppm)} ppm")
    print(f"    희석배수 : {PREP.dilution_factor:.0f}배 "
          f"(원액 {PREP.drink_aliquot_mL:.0f} mL → {PREP.final_volume_mL:.0f} mL)")
    print()
    from .simulate import TRUE_DRINK_PPM

    for drink_key, d in DRINKS.items():
        print(f"  [{d.name_ko}] (문헌 일반값 기준 가정)")
        expected = {
            k: TRUE_DRINK_PPM[drink_key][k] / PREP.dilution_factor
            for k in COMPOUND_ORDER
        }
        for line in an.check_spike_design(expected):
            print(f"    {line}")
        print()

    print("[3] 표준액 조제량")
    for spike in PREP.spike_levels_ppm:
        v = PREP.spike_volume_mL(spike)
        print(f"    +{spike:>3.0f} ppm  →  {PREP.stock_ppm:.0f} ppm 표준액 {v:.2f} mL")
    print()
    n_cmp = len(COMPOUND_ORDER)
    n_vials = len(DRINKS) * len(PREP.spike_levels_ppm) + len(PREP.spike_levels_ppm)
    print("[4] 놓치기 쉬운 것")
    print(f"    · 표준액을 성분별로 따로 만들면 바이알 수가 {n_cmp}배가 됩니다.")
    print(f"      {n_cmp}성분 혼합 표준액 1종으로 만들면 한 세트로 끝납니다.")
    print(f"      (혼합 표준액 기준 총 바이알 {n_vials}개 = 검량선 "
          f"{len(PREP.spike_levels_ppm)} + 음료 {len(DRINKS)}종 x "
          f"{len(PREP.spike_levels_ppm)})")
    print("    · 모든 바이알의 음료 원액 양(1 mL)은 정확히 같아야 합니다.")
    print("      이게 어긋나면 x절편 역산 전제가 깨집니다.")
    print("    · 탄산 제거 전후로 부피가 변하므로, 탈기 후에 분취하세요.")
    print("    · 시린지 필터는 첫 1 mL 를 버리고 받으세요(필터 흡착 손실).")
    return 0


def cmd_simulate(args) -> int:
    from .simulate import SimulationSettings, write_simulated_csv

    st = SimulationSettings(
        seed=args.seed,
        include_matrix=not args.no_matrix,
        include_additive=not args.no_additive,
    )
    path = write_simulated_csv(args.output, st)
    print(f"모의 데이터를 만들었습니다: {path}")
    print("⚠ 실측값이 아닙니다. 코드 검증 및 예상 결과 확인용입니다.")
    print(f"다음: python -m hplc analyze {path}")
    return 0


def cmd_chromatogram(args) -> int:
    from .simulate import SimulationSettings, simulate_chromatogram

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plots.setup_style()
    st = SimulationSettings(seed=args.seed)

    made = []
    t, sig, peaks = simulate_chromatogram(None, spike_ppm=40.0, settings=st)
    made.append(
        plots.plot_chromatogram(
            t, sig, peaks,
            plots.L("혼합 표준액 40 ppm (증류수 바탕) — 모의",
                    "Mixed standard 40 ppm (aqueous) - simulated"),
            outdir, "fig_chromatogram_standard.png",
        )
    )
    for drink_key, d in DRINKS.items():
        t, sig, peaks = simulate_chromatogram(drink_key, spike_ppm=0.0, settings=st)
        made.append(
            plots.plot_chromatogram(
                t, sig, peaks,
                plots.L(f"{d.name_ko} 무첨가 시료 — 모의",
                        f"{drink_key} unspiked - simulated"),
                outdir, f"fig_chromatogram_{drink_key}.png",
            )
        )
    for p in made:
        print(f"  {p}")
    print("⚠ 모의 크로마토그램입니다. 실측 크로마토그램 대신 쓰지 마세요.")
    return 0


def cmd_analyze(args) -> int:
    src = Path(args.csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    simulated = dataio.is_simulated(src)
    try:
        df = dataio.load_peak_areas(src)
    except dataio.DataValidationError as e:
        print(f"입력 오류: {e}", file=sys.stderr)
        return 1

    if simulated:
        print("=" * 68)
        print("⚠ 모의(SIMULATED) 데이터를 분석합니다. 결과는 실측값이 아닙니다.")
        print("=" * 68)

    warns = dataio.check_retention_times(df)

    try:
        cals = an.build_calibrations(df)
        results = an.run_standard_addition(df, cals)
    except ValueError as e:
        print(f"분석 오류: {e}", file=sys.stderr)
        return 1

    if not results:
        print("분석할 시료 데이터가 없습니다.", file=sys.stderr)
        return 1

    figures: list[Path] = []
    if not args.no_plots:
        plots.setup_style()
        cal_fig = plots.plot_calibrations(cals, outdir)
        if cal_fig is not None:      # 검량선을 계수로만 넣은 경우 그릴 점이 없다
            figures.append(cal_fig)
        figures.extend(plots.plot_standard_addition(results, outdir))
        conv_fig = plots.plot_conversion_constants(results, outdir)
        if conv_fig is not None:
            figures.append(conv_fig)
        figures.append(plots.plot_residuals(results, outdir))

    warns.extend(an.calibration_missing_note(cals))

    md = report.build_report(
        cals, results,
        simulated=simulated, source=str(src),
        figures=figures, extra_warnings=warns,
    )
    prefix = "SIMULATED_" if simulated else ""
    md_path = outdir / f"{prefix}report.md"
    md_path.write_text(md, encoding="utf-8")

    res_df = report.results_dataframe(results)
    cal_df = report.calibration_dataframe(cals)
    res_path = outdir / f"{prefix}results_summary.csv"
    cal_path = outdir / f"{prefix}calibration_summary.csv"
    res_df.to_csv(res_path, index=False, encoding="utf-8-sig")
    cal_df.to_csv(cal_path, index=False, encoding="utf-8-sig")

    # 콘솔 요약
    print()
    print("성분별 역산 결과 (음료 원액 기준)")
    print("-" * 68)
    for r in results:
        import math as _math
        f_txt = (f"f={r.conversion_constant:>6.3f}"
                 if _math.isfinite(r.conversion_constant) else "f=     -")
        print(
            f"  {r.drink.name_ko:<16} {r.compound.name_ko:<10} "
            f"{r.c_drink_ppm:>8.1f} ppm   "
            f"{r.mg_per_serving:>7.1f} mg/캔   "
            f"{f_txt}   R²={r.fit.r2:.4f}"
        )
    print("-" * 68)
    print()
    print("생성된 파일")
    print(f"  보고서 : {md_path}")
    print(f"  결과표 : {res_path}")
    print(f"  검량선 : {cal_path}")
    for p in figures:
        print(f"  그림   : {p}")

    all_warnings = list(warns) + [w for r in results for w in r.warnings()]
    if all_warnings:
        print()
        print(f"점검 필요 항목 {len(set(all_warnings))}건 — 보고서 하단을 확인하세요.")
    return 0


def cmd_demo(args) -> int:
    from .simulate import SimulationSettings, write_simulated_csv

    DEFAULT_DATA.mkdir(parents=True, exist_ok=True)
    csv_path = DEFAULT_DATA / "SIMULATED_peak_areas.csv"
    write_simulated_csv(csv_path, SimulationSettings(seed=args.seed))
    print(f"[1/3] 모의 데이터 생성: {csv_path}")

    ns = argparse.Namespace(outdir=str(DEFAULT_OUT), seed=args.seed)
    print("[2/3] 모의 크로마토그램 생성")
    cmd_chromatogram(ns)

    print("[3/3] 분석 실행")
    ns2 = argparse.Namespace(csv=str(csv_path), outdir=str(DEFAULT_OUT), no_plots=False)
    return cmd_analyze(ns2)


# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hplc",
        description="제로 에너지 드링크 HPLC 동시분석 - 표준물 첨가법 정량 도구",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = p.add_subparsers(dest="command", required=True)

    t = sub.add_parser("template", help="빈 입력 CSV 템플릿 생성")
    t.add_argument("-o", "--output", default=str(DEFAULT_DATA / "peak_areas_template.csv"))
    t.set_defaults(func=cmd_template)

    d = sub.add_parser("design-check", help="실험 전 설계 점검")
    d.set_defaults(func=cmd_design_check)

    s = sub.add_parser("simulate", help="모의 데이터 CSV 생성 (검증용)")
    s.add_argument("-o", "--output", default=str(DEFAULT_DATA / "SIMULATED_peak_areas.csv"))
    s.add_argument("--seed", type=int, default=20260812)
    s.add_argument("--no-matrix", action="store_true", help="곱셈형 매트릭스 효과 끄기")
    s.add_argument("--no-additive", action="store_true", help="덧셈형 공용리 간섭 끄기")
    s.set_defaults(func=cmd_simulate)

    c = sub.add_parser("chromatogram", help="모의 크로마토그램 그림 생성")
    c.add_argument("-o", "--outdir", default=str(DEFAULT_OUT))
    c.add_argument("--seed", type=int, default=20260812)
    c.set_defaults(func=cmd_chromatogram)

    a = sub.add_parser("analyze", help="CSV 분석 -> 보고서 + 그림")
    a.add_argument("csv")
    a.add_argument("-o", "--outdir", default=str(DEFAULT_OUT))
    a.add_argument("--no-plots", action="store_true")
    a.set_defaults(func=cmd_analyze)

    m = sub.add_parser("demo", help="모의 생성부터 분석까지 한 번에")
    m.add_argument("--seed", type=int, default=20260812)
    m.set_defaults(func=cmd_demo)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
