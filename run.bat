@echo off
chcp 65001 >nul
title 동심 이중관 층류 유동 시뮬레이터 (Concentric Annular Flow)
cd /d "%~dp0"

echo ============================================================
echo   🌊 동심 이중관 층류 유동 시뮬레이터 자동 실행기
echo ============================================================
echo.

:: 1. 시스템에 설치된 Python 탐색
set PYTHON_CMD=

python --version >nul 2>&1
if %errorlevel% == 0 (
    set PYTHON_CMD=python
    goto FOUND_PYTHON
)

py -3 --version >nul 2>&1
if %errorlevel% == 0 (
    set PYTHON_CMD=py -3
    goto FOUND_PYTHON
)

py --version >nul 2>&1
if %errorlevel% == 0 (
    set PYTHON_CMD=py
    goto FOUND_PYTHON
)

:: 기본 사용자 경로 탐색 (이전 환경 호환)
if exist "%LOCALAPPDATA%\Programs\Python\Python312\python.exe" (
    set "PYTHON_CMD=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
    goto FOUND_PYTHON
)

if exist "C:\Users\USER\AppData\Local\Programs\Python\Python312\python.exe" (
    set "PYTHON_CMD=C:\Users\USER\AppData\Local\Programs\Python\Python312\python.exe"
    goto FOUND_PYTHON
)

echo [❌ 오류] 컴퓨터에 Python이 설치되어 있지 않거나 PATH 환경변수에 등록되지 않았습니다.
echo.
echo 💡 해결 방법:
echo    1) https://www.python.org/downloads/ 에서 Python 최신 버전을 설치해 주세요.
echo    2) 설치 시 반드시 화면 하단의 "Add python.exe to PATH" 옵션을 체크해야 합니다.
echo.
pause
exit /b 1

:FOUND_PYTHON
echo [1/3] Python 환경 감지 완료: %PYTHON_CMD%

:: 2. 필수 라이브러리 검사 및 자동 설치
echo [2/3] 필수 라이브러리 (Streamlit, NumPy, Matplotlib, Plotly, SciPy) 확인 및 설치 중...
%PYTHON_CMD% -m pip install -r requirements.txt --quiet --disable-pip-version-check

:: 3. 웹 시뮬레이터 구동
echo [3/3] 웹 시뮬레이션 서버를 시작합니다...
echo.
echo ============================================================
echo  🌐 웹 브라우저가 잠시 후 자동으로 열립니다.
echo.
echo  * 로컬 접속 주소:     http://localhost:8501
echo  * 같은 Wi-Fi 접속:    터미널에 표시되는 Network URL (예: http://192.168.x.x:8501)
echo  * 프로그램 종료 방법:  이 창에서 Ctrl + C
echo ============================================================
echo.

%PYTHON_CMD% -m streamlit run app.py

pause
