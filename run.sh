#!/usr/bin/env bash
# 동심 이중관 층류 유동 시뮬레이터 (macOS / Linux 전용 실행 스크립트)
cd "$(dirname "$0")"

echo "============================================================"
echo " 🌊 동심 이중관 층류 유동 시뮬레이터 (macOS / Linux)"
echo "============================================================"

# Python3 명령어 탐색
if command -v python3 &>/dev/null; then
    PY_CMD=python3
elif command -v python &>/dev/null; then
    PY_CMD=python
else
    echo "[오류] Python이 설치되어 있지 않습니다."
    exit 1
fi

echo "[1/2] 필수 라이브러리 설치 확인 중..."
$PY_CMD -m pip install -r requirements.txt --quiet

echo "[2/2] 시뮬레이션 서버 실행 중..."
echo " * 로컬 주소: http://localhost:8501"
echo " * 종료 방법: Ctrl + C"
echo "============================================================"

$PY_CMD -m streamlit run app.py
