"""colab.ipynb 를 저장소 소스로부터 생성한다.

노트북이 저장소 내려받기 없이 혼자 돌아가야 하므로, 프로젝트 파일 전체를
zip 으로 묶어 base64 로 노트북 안에 심는다. 노트북 두 번째 셀이 그것을
/content/hplc_project 에 풀어 놓으면 그 뒤로는 평소처럼 쓰면 된다.

    python tools/build_colab.py

코드를 고쳤으면 이 스크립트를 다시 돌려서 colab.ipynb 를 갱신할 것.
(그렇게 하지 않으면 노트북만 옛 코드를 담은 채로 남는다.)
"""

from __future__ import annotations

import base64
import io
import json
import textwrap
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "colab.ipynb"

# 노트북에 심을 파일들
INCLUDE = [
    "hplc/__init__.py",
    "hplc/__main__.py",
    "hplc/analysis.py",
    "hplc/cli.py",
    "hplc/config.py",
    "hplc/dataio.py",
    "hplc/plots.py",
    "hplc/report.py",
    "hplc/simulate.py",
    "hplc/stats.py",
    "tests/test_analysis.py",
    "data/measured_peak_areas.csv",
    "requirements.txt",
]

PROJECT_DIR = "/content/hplc_project"


def build_payload() -> tuple[str, int]:
    """프로젝트 파일을 zip -> base64 로 만든다."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for rel in INCLUDE:
            src = ROOT / rel
            if not src.exists():
                raise SystemExit(f"파일이 없습니다: {rel}")
            zf.writestr(rel, src.read_bytes())
    raw = buf.getvalue()
    return base64.b64encode(raw).decode("ascii"), len(raw)


def b64_cell_lines(b64: str) -> list[str]:
    """base64 문자열을 읽기 좋은 폭으로 잘라 파이썬 리터럴로 만든다."""
    chunks = textwrap.wrap(b64, 100)
    lines = ["PROJECT_ZIP_B64 = (\n"]
    lines += [f'    "{c}"\n' for c in chunks]
    lines.append(")\n")
    return lines


def md(*text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": list(text)}


def code(*text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": list(text),
    }


GUARD = (
    "import os\n"
    f"os.chdir('{PROJECT_DIR}')   # 셀 순서가 어긋나도 되도록 매번 이동\n"
    "\n"
)


def build_notebook() -> dict:
    b64, raw_size = build_payload()

    cells = [
        md(
            "# 제로 에너지 드링크 HPLC 분석\n",
            "\n",
            "카페인 · 소듐벤조에이트를 표준물 첨가법으로 정량합니다.\n",
            "\n",
            "**메뉴에서 `런타임 → 모두 실행` 을 누르면 끝까지 돌아갑니다.**\n",
            "따로 내려받거나 업로드할 파일은 없습니다. 코드와 실측 데이터가\n",
            "이 노트북 안에 들어 있습니다.\n",
            "\n",
            "---\n",
            "\n",
            "### ⚠ 모의 데이터에 관하여\n",
            "\n",
            "맨 아래 선택 항목에 **모의 데이터 생성기**가 있습니다. 코드가 제대로\n",
            "도는지 확인하는 용도입니다. 거기서 나온 숫자를 보고서에 실측값으로\n",
            "적으면 **데이터 조작**입니다. 모의 데이터로 돌린 결과물에는 `SIMULATED`\n",
            "표시가 자동으로 박히니 지우지 마세요.\n",
        ),
        md(
            "## 1단계 — 준비물 설치\n",
            "\n",
            "분석 라이브러리와 한글 폰트를 깝니다. 1~2분 걸립니다.\n",
            "\n",
            "폰트를 까는 이유는 matplotlib이 한글 폰트 없이 그래프를 그리면 글자가\n",
            "네모(□□□)로 나오기 때문입니다.\n",
        ),
        code(
            "!pip install -q numpy pandas matplotlib scipy\n",
            "!apt-get install -y -qq fonts-nanum > /dev/null 2>&1\n",
            "!fc-cache -f > /dev/null 2>&1\n",
            "\n",
            "import matplotlib.font_manager as fm\n",
            "fm._load_fontmanager(try_read_cache=False)\n",
            "print('한글 폰트:', 'NanumGothic' in {f.name for f in fm.fontManager.ttflist})\n",
            "print('설치 완료')\n",
        ),
        md(
            "## 2단계 — 코드 풀어놓기\n",
            "\n",
            f"분석 코드와 실측 데이터를 `{PROJECT_DIR}` 에 풉니다.\n",
            "인터넷에서 받아오지 않고 이 노트북 안에 심어 둔 것을 씁니다.\n",
            "\n",
            "푼 뒤에는 왼쪽 폴더 아이콘에서 파일을 직접 열어 고칠 수 있습니다.\n",
        ),
        code(
            "import base64, io, zipfile, os, sys, shutil\n",
            "\n",
            *b64_cell_lines(b64),
            "\n",
            f"PROJECT = '{PROJECT_DIR}'\n",
            "shutil.rmtree(PROJECT, ignore_errors=True)\n",
            "os.makedirs(PROJECT, exist_ok=True)\n",
            "with zipfile.ZipFile(io.BytesIO(base64.b64decode(PROJECT_ZIP_B64))) as zf:\n",
            "    zf.extractall(PROJECT)\n",
            "os.chdir(PROJECT)\n",
            "if PROJECT not in sys.path:\n",
            "    sys.path.insert(0, PROJECT)\n",
            "\n",
            "for root, _, names in os.walk(PROJECT):\n",
            "    for n in sorted(names):\n",
            "        print(os.path.relpath(os.path.join(root, n), PROJECT))\n",
        ),
        md(
            "## 3단계 — 실험 설계 점검\n",
            "\n",
            "머무름 시간, 머무름 계수 k', 첨가 농도 적정성, 표준액 조제량을 출력합니다.\n",
            "파일은 만들지 않고 화면에만 나옵니다.\n",
        ),
        code(GUARD, "!python -m hplc design-check\n"),
        md(
            "## 4단계 — 실측 데이터 분석\n",
            "\n",
            "실험에서 나온 피크 면적 96행(음료 3종 × 성분 2종 × 첨가 4수준 × 주입 4회)을\n",
            "분석합니다. 표준물 첨가법으로 x절편을 역산해 실제 농도를 구합니다.\n",
            "\n",
            "> 증류수 바탕 검량선 데이터가 아직 없어서 **변환 상수 칸은 비어 나옵니다.**\n",
            "> 농도 역산은 검량선을 쓰지 않으므로 결과 자체는 그대로 유효합니다.\n",
            "> 검량선을 구하셨으면 6단계 아래의 안내를 보세요.\n",
        ),
        code(GUARD, "!python -m hplc analyze data/measured_peak_areas.csv\n"),
        md(
            "## 5단계 — 그림 보기\n",
            "\n",
            "표준물 첨가법 그래프와 잔차 그래프를 노트북 안에서 봅니다.\n",
            "그래프를 우클릭하면 이미지로 저장할 수 있습니다.\n",
        ),
        code(
            GUARD,
            "from IPython.display import Image, display, Markdown\n",
            "import glob\n",
            "\n",
            "figs = sorted(glob.glob('output/*.png'))\n",
            "if not figs:\n",
            "    print('그림이 없습니다. 4단계를 먼저 실행하세요.')\n",
            "for f in figs:\n",
            "    display(Markdown(f'### {os.path.basename(f)}'))\n",
            "    display(Image(filename=f))\n",
        ),
        md(
            "## 6단계 — 보고서 읽기\n",
            "\n",
            "보고서 7장에 그대로 옮길 수 있는 표들입니다. 자동 점검 결과도 맨 아래\n",
            "붙어 있으니 꼭 확인하세요.\n",
        ),
        code(
            GUARD,
            "from IPython.display import Markdown, display\n",
            "import glob\n",
            "\n",
            "reports = sorted(glob.glob('output/*report.md'))\n",
            "if reports:\n",
            "    display(Markdown(open(reports[0], encoding='utf-8').read()))\n",
            "else:\n",
            "    print('보고서가 없습니다. 4단계를 먼저 실행하세요.')\n",
        ),
        md(
            "## 7단계 — 결과 내려받기\n",
            "\n",
            "보고서·표·그림을 ZIP 하나로 묶어 저장합니다.\n",
        ),
        code(
            GUARD,
            "import shutil\n",
            "from google.colab import files\n",
            "\n",
            "if not os.path.isdir('output'):\n",
            "    print('결과가 없습니다. 4단계를 먼저 실행하세요.')\n",
            "else:\n",
            "    shutil.make_archive('hplc_results', 'zip', 'output')\n",
            "    files.download('hplc_results.zip')\n",
        ),
        md(
            "---\n",
            "\n",
            "# 여기부터는 선택 항목\n",
            "\n",
            "`런타임 → 모두 실행` 은 위 7단계까지로 충분합니다. 아래는 필요할 때만\n",
            "개별로 실행하세요.\n",
            "\n",
            "## 값 고치기 — 라벨 표시량, 검량선, 실험 조건\n",
            "\n",
            "왼쪽 **폴더 아이콘 → hplc → config.py** 를 더블클릭하면 편집기가 열립니다.\n",
            "고친 뒤 4단계부터 다시 실행하면 반영됩니다.\n",
            "\n",
            "| 고칠 곳 | 무엇 |\n",
            "|---|---|\n",
            "| `DRINKS` | 캔의 실제 카페인 표기량(mg)과 용량(mL) ← **지금 자리표시자입니다** |\n",
            "| `EXTERNAL_CALIBRATION` | 검량선 계수 `{\"caffeine\": (기울기, y절편), ...}` |\n",
            "| `HPLCConditions` | 이동상 비율, 유량, 오븐 온도, 검출 파장 |\n",
            "| `expected_rt_min` | 실측 머무름 시간 |\n",
            "| `PrepConditions` | 희석배수, 첨가 농도 |\n",
            "\n",
            "> Colab에서 고친 내용은 런타임이 끊기면 사라집니다. 계속 쓸 값이면\n",
            "> GitHub의 `config.py` 도 같이 고쳐 두세요.\n",
        ),
        md(
            "## (선택) 내 CSV 올려서 분석하기\n",
            "\n",
            "검량선 행을 추가했거나 재측정한 데이터가 있으면 여기에 올리세요.\n",
            "실행하면 파일 선택 버튼이 나옵니다.\n",
        ),
        code(
            GUARD,
            "import shutil\n",
            "from google.colab import files\n",
            "\n",
            "shutil.rmtree('output', ignore_errors=True)\n",
            "print('CSV 파일을 선택하세요.')\n",
            "uploaded = files.upload()\n",
            "csv_name = next(n for n in uploaded if n.lower().endswith('.csv'))\n",
            "\n",
            "# Colab 은 같은 이름이 있으면 'xxx (1).csv' 로 저장하므로,\n",
            "# 올린 내용을 항상 같은 이름으로 직접 써 준다.\n",
            "MEASURED = 'uploaded.csv'\n",
            "with open(MEASURED, 'wb') as fh:\n",
            "    fh.write(uploaded[csv_name])\n",
            "print(f'{csv_name} -> {MEASURED} ({len(uploaded[csv_name]):,} bytes)')\n",
            "\n",
            "!python -m hplc analyze uploaded.csv\n",
        ),
        md(
            "## (선택) 빈 입력표 받기\n",
            "\n",
            "새로 실험할 때 쓸 빈 CSV입니다. 구글 스프레드시트로 열어\n",
            "`peak_area` 열을 채우세요. 피크가 없으면 `0` 을 적습니다.\n",
        ),
        code(
            GUARD,
            "from google.colab import files\n",
            "\n",
            "!python -m hplc template\n",
            "files.download('data/peak_areas_template.csv')\n",
        ),
        md(
            "## (선택) 모의 데이터로 코드 검증\n",
            "\n",
            "정답을 아는 가짜 데이터를 넣고 그 정답이 되돌아오는지 확인합니다.\n",
            "**여기서 나오는 숫자는 실측값이 아닙니다.**\n",
        ),
        code(GUARD, "!python -m pytest tests/ -q\n", "!python -m hplc demo\n"),
    ]

    return {
        "cells": cells,
        "metadata": {
            "colab": {"provenance": [], "toc_visible": True},
            "kernelspec": {"display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 0,
    }, raw_size


def main() -> None:
    nb, raw_size = build_notebook()
    with OUT.open("w", encoding="utf-8") as fh:
        json.dump(nb, fh, ensure_ascii=False, indent=1)
    size = OUT.stat().st_size
    print(f"{OUT.name} 생성 완료")
    print(f"  포함 파일 {len(INCLUDE)}개, 압축 전 {raw_size:,} bytes")
    print(f"  노트북 크기 {size:,} bytes, 셀 {len(nb['cells'])}개")


if __name__ == "__main__":
    main()
