@echo off
chcp 65001 >nul
title 동심 이중관 층류 유동 시뮬레이터 (가상환경 완전 격리 실행기)
cd /d "%~dp0"

echo ============================================================
echo   🌊 동심 이중관 층류 유동 시뮬레이터 (가상환경 venv 실행기)
echo ============================================================
echo.

:: 1. Python 탐색
set PYTHON_CMD=
python --version >nul 2>&1 && set PYTHON_CMD=python
if "%PYTHON_CMD%"=="" py -3 --version >nul 2>&1 && set PYTHON_CMD=py -3
if "%PYTHON_CMD%"=="" py --version >nul 2>&1 && set PYTHON_CMD=py
if "%PYTHON_CMD%"=="" (
    echo [❌ 오류] Python이 설치되어 있지 않습니다. python.org에서 설치해 주세요.
    pause
    exit /b 1
)

:: 2. 가상환경 생성 여부 확인
if not exist "venv\Scripts\activate.bat" (
    echo [1/3] 독립 가상환경(venv)을 처음 생성합니다...
    %PYTHON_CMD% -m venv venv
    if errorlevel 1 (
        echo [❌ 오류] 가상환경 생성에 실패했습니다.
        pause
        exit /b 1
    )
    echo [2/3] 가상환경에 필수 패키지를 설치합니다 (최초 1회, 1~2분 소요)...
    call venv\Scripts\activate.bat
    python -m pip install --upgrade pip --quiet
    python -m pip install -r requirements.txt --quiet
) else (
    echo [1/2] 기존 가상환경(venv) 활성화...
    call venv\Scripts\activate.bat
)

:: 3. 실행
echo [2/2] 시뮬레이션 서버 구동 중...
echo.
echo ============================================================
echo  🌐 웹 브라우저가 잠시 후 자동으로 열립니다.
echo  * 로컬 접속:   http://localhost:8501
echo  * 종료 방법:   이 창에서 Ctrl + C
echo ============================================================
echo.

streamlit run app.py

pause
