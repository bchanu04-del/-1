# 동심 이중관 층류 유동 시뮬레이터 범용 PowerShell 실행 스크립트
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $scriptDir

# Python 실행 파일 탐색
$pythonCmd = $null
if (Get-Command python.exe -ErrorAction SilentlyContinue) {
    $pythonCmd = "python"
} elseif (Get-Command py.exe -ErrorAction SilentlyContinue) {
    $pythonCmd = "py"
} elseif (Test-Path "$env:LOCALAPPDATA\Programs\Python\Python312\python.exe") {
    $pythonCmd = "$env:LOCALAPPDATA\Programs\Python\Python312\python.exe"
}

if (-not $pythonCmd) {
    Write-Host "[❌ 오류] 컴퓨터에 Python이 설치되어 있지 않거나 PATH 환경변수에 등록되지 않았습니다." -ForegroundColor Red
    Write-Host "https://www.python.org/downloads/ 에서 Python을 설치해 주세요. (설치 시 Add to PATH 필수 체크)" -ForegroundColor Yellow
    pause
    exit 1
}

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host " 🌊 동심 이중관 층류 유동 시뮬레이터 (Concentric Annular Flow)" -ForegroundColor Green
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "[1/2] 필수 라이브러리 검사 및 설치..." -ForegroundColor Gray
& $pythonCmd -m pip install -r requirements.txt --quiet --disable-pip-version-check

Write-Host "[2/2] 시뮬레이션 서버 구동 중..." -ForegroundColor Gray
Write-Host " * 로컬 접속:   http://localhost:8501" -ForegroundColor Yellow
Write-Host " * 종료 방법:   터미널에서 Ctrl + C`n" -ForegroundColor DarkGray

& $pythonCmd -m streamlit run app.py
