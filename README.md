# 동심 이중관 내 층류 유동 시뮬레이션 (Concentric Annular Laminar Flow Simulation)

화학공학 유체역학 이론에 기반하여 동심 이중관(Concentric Annulus) 사이에서 발생하는 정상 상태(Steady-state)의 완전 발달(Fully-developed) 층류 유동을 수치 계산하고, 유속 및 전단응력 분포를 웹 UI로 시각화하는 시뮬레이션 프로그램입니다.

---

## 1. 주요 기능 및 특징

- **엄밀한 해석해(Analytical Solution) 기반 계산**: 원통 좌표계 Navier–Stokes 방정식을 엄밀하게 푼 유속 분포 $u(r)$ 및 전단응력 분포 $\tau_{rz}(r)$ 산출.
- **실시간 반응형 웹 인터페이스 (Streamlit)**: 파이프 크기, 압력차, 유체 물성 변경 시 지표 및 그래프 실시간 자동 갱신.
- **다양한 시각화 지원**:
  - Matplotlib 고해상도 출판용 차트 탭
  - Plotly 인터랙티브 차트 탭 (마우스 호버 수치 확인, 줌, 패닝 지원)
  - 2차원 환형 단면(Cross-section) 등고선도
- **수치적 이론 검증 패널**:
  - 벽면 No-slip 경계조건 ($u(R_i)=0, u(R_o)=0$) 만족 검증 ($< 10^{-12}$ m/s)
  - 최대속도 위치($r_{max}$)에서의 속도 구배 $\frac{du}{dr} = 0$ 및 전단응력 $\tau_{rz} = 0$ 검증
  - 엄밀 해석해 체적 유량($Q_{analytical}$)과 사다리꼴 수치 적분 유량($Q_{numerical}$)의 일치성 검증 (상대 오차 < 0.01%)
  - 레이놀즈 수($Re$) 계산 및 층류/천이/난류 판정 가이드 제공
- **벽면 전단응력 부호 및 방향 해석**: 내부벽($\tau > 0$)과 외부벽($\tau < 0$)의 물리적 마찰 방향 설명.
- **어디서나 실행 가능한 포터블 구조**: USB나 압축 파일로 다른 컴퓨터에 옮겨도 더블 클릭 한 번으로 자동 라이브러리 설치 및 실행 지원.

---

## 2. 이론적 배경

### A. 원통 좌표계 Navier–Stokes 지배 방정식
정상 상태, 비압축성 뉴턴 유체, 축대칭 유동($\partial/\partial\theta = 0$), 축방향 완전 발달 유동($v_r = v_\theta = 0, \partial u/\partial z = 0$) 가정:

$$0 = -\frac{dp}{dz} + \mu \frac{1}{r} \frac{d}{dr}\left(r \frac{du}{dr}\right)$$

여기서 압력 강하는 $\Delta P = P_{in} - P_{out} > 0$ 이며, $-\frac{dp}{dz} = \frac{\Delta P}{L}$ 입니다.

### B. 유속 분포 (Velocity Profile)
이중 No-slip 경계조건 $u(R_i) = 0$, $u(R_o) = 0$을 적용하여 적분:

$$u(r) = \frac{\Delta P}{4\mu L} \left[ R_o^2 - r^2 - \frac{R_o^2 - R_i^2}{\ln(R_o / R_i)} \ln\left(\frac{R_o}{r}\right) \right]$$

### C. 최대 유속 위치 및 최대 유속
$\frac{du}{dr} = 0$ 인 극값 조건으로부터:

$$r_{max} = \sqrt{\frac{R_o^2 - R_i^2}{2 \ln(R_o / R_i)}}, \quad u_{max} = u(r_{max})$$

> **물리적 고찰**: 원통 좌표계 곡률로 인해 $r_{max}$는 두 관의 산술 평균 $\frac{R_i + R_o}{2}$보다 항상 안쪽(내부관 쪽)으로 치우칩니다.

### D. 전단응력 분포 (Shear Stress Distribution)
$$\tau_{rz}(r) = \mu \frac{du}{dr} = \frac{\Delta P}{4L} \left[ -2r + \frac{R_o^2 - R_i^2}{r \ln(R_o / R_i)} \right]$$

- 내부관 벽면 ($r = R_i$): $\frac{du}{dr} > 0 \implies \tau_{rz}(R_i) > 0$
- 외부관 벽면 ($r = R_o$): $\frac{du}{dr} < 0 \implies \tau_{rz}(R_o) < 0$
- $r = r_{max}$ 지점: $\frac{du}{dr} = 0 \implies \tau_{rz}(r_{max}) = 0$

### E. 유량, 평균 유속 및 레이놀즈 수
- **환형 단면적**: $A = \pi(R_o^2 - R_i^2)$
- **수력 직경**: $D_h = 2(R_o - R_i)$
- **체적 유량 (해석해)**:
  $$Q = 2\pi \int_{R_i}^{R_o} u(r)r\,dr = \frac{\pi \Delta P}{8\mu L} \left[ (R_o^4 - Ri^4) - \frac{(R_o^2 - R_i^2)^2}{\ln(R_o / R_i)} \right]$$
- **평균 유속**: $u_{avg} = \frac{Q}{A}$
- **레이놀즈 수**: $Re = \frac{\rho u_{avg} D_h}{\mu}$

---

## 3. 프로젝트 파일 구성

```
concentric_annular_flow/
├── app.py               # Streamlit 기반 인터랙티브 웹 UI 및 시각화 대시보드
├── flow_model.py        # 유속, 전단응력, 최대속도, 유량, Re 계산 및 검증 모델 클래스
├── test_flow_model.py   # pytest 기반 단위 테스트 (경계조건, 극값, 오차 검증)
├── requirements.txt     # 의존성 라이브러리 목록 (Streamlit, NumPy, Matplotlib, Plotly, SciPy, pytest)
├── run.bat              # [Windows] 더블 클릭 원클릭 자동 설치 및 실행기
├── run_with_venv.bat    # [Windows] 독립 가상환경(venv) 자동 생성 및 실행기
├── run.ps1              # [PowerShell] 범용 실행 스크립트
├── run.sh               # [macOS / Linux] 실행 스크립트
├── .streamlit/          # Streamlit 바인딩 설정 파일
└── README.md            # 프로그램 설명서 및 사용 가이드
```

---

## 4. 다른 컴퓨터/노트북에서 실행하는 방법

### 방법 A. 프로젝트 폴더를 다른 컴퓨터로 복사하여 실행 (가장 추천)
1. 생성된 **`concentric_annular_flow_portable.zip`** 파일을 USB, 카카오톡, 이메일, Google Drive 등으로 다른 노트북에 전송합니다.
2. 다른 노트북에서 압축을 풉니다.
3. 압축 푼 폴더에서 **`run.bat`** 파일을 **더블 클릭**합니다.
   - 필요한 라이브러리(NumPy, Matplotlib 등)가 자동으로 감지되어 자동 설치된 후 브라우저가 열립니다.
   - (참고: 다른 컴퓨터에 Python 3.10 이상이 설치되어 있어야 합니다.)

### 방법 B. 같은 Wi-Fi에 있는 다른 컴퓨터/스마트폰에서 접속하기
현재 노트북에서 프로그램을 켜둔 상태라면, 동일한 공유기(Wi-Fi)에 연결된 다른 노트북이나 태블릿, 스마트폰에서 프로그램 설치 없이 웹 브라우저로 바로 조작할 수 있습니다.
1. 현재 노트북에서 `run.bat`을 실행합니다.
2. 콘솔 창에 표시되는 **Network URL**(예: `http://192.168.0.15:8501`)을 확인합니다.
3. 다른 노트북이나 스마트폰의 브라우저 주소창에 해당 주소를 입력하면 바로 접속됩니다.

### 방법 C. 인터넷에 무료 배포하여 어디서나 접속하기 (Streamlit Community Cloud)
GitHub 계정이 있다면 단 2분 만에 인터넷 무료 웹사이트로 배포할 수 있습니다:
1. GitHub에 새 저장소(Public)를 만들고 이 폴더의 파일들을 업로드합니다.
2. [share.streamlit.io](https://share.streamlit.io)에 접속하여 GitHub 저장소와 `app.py`를 선택합니다.
3. 무료 고유 URL(예: `https://my-flow-sim.streamlit.app`)이 발급되어 전 세계 어떤 기기에서든 링크만 누르면 바로 사용할 수 있습니다.

---

## 5. 입력 변수 설명

| 변수명 | 기호 | 기본 단위 | 유효 범위 | 설명 |
|---|---|---|---|---|
| 내부관 반지름 | $R_i$ | mm (내부 변환 m) | $R_i > 0$ | 내부관 바깥쪽 반지름 |
| 외부관 내벽 반지름 | $R_o$ | mm (내부 변환 m) | $R_o > R_i$ | 외부관 안쪽 반지름 |
| 관 길이 | $L$ | m | $L > 0$ | 배관의 축방향 전체 길이 |
| 유체 점도 | $\mu$ | Pa·s | $\mu > 0$ | 유체의 역학적 점성계수 (물: 약 0.001 Pa·s) |
| 유체 밀도 | $\rho$ | kg/m³ | $\rho > 0$ | 유체 밀도 (물: 약 998 kg/m³) |
| 압력 강하 | $\Delta P$ | Pa | $\Delta P > 0$ | 입구 압력과 출구 압력의 차이 ($P_{in} - P_{out}$) |
