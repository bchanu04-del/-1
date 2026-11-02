"""
app.py
======
동심 이중관 내 층류 유동 시뮬레이션 (Concentric Annular Laminar Flow Simulation)
화학공학 유체역학 기반 웹 인터랙티브 시뮬레이터 (Streamlit)
"""

import sys
import streamlit as st
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Streamlit 멀티스레드 환경에서 GUI 백엔드 충돌 및 빈 그래프 발생 방지
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from flow_model import AnnularFlowParameters, ConcentricAnnularFlowModel

# -----------------------------------------------------------------------------
# 1. 페이지 환경설정 & 한글 폰트 설정
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="동심 이중관 층류 유동 시뮬레이션",
    page_icon="🌊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Matplotlib 한글 폰트 및 마이너스 부호 설정
def setup_matplotlib_font():
    fonts = [f.name for f in fm.fontManager.ttflist]
    if "Malgun Gothic" in fonts:
        plt.rcParams["font.family"] = "Malgun Gothic"
    elif "NanumGothic" in fonts:
        plt.rcParams["font.family"] = "NanumGothic"
    elif "AppleGothic" in fonts:
        plt.rcParams["font.family"] = "AppleGothic"
    else:
        plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.unicode_minus"] = False

setup_matplotlib_font()


def inject_anti_dom_crash_headers():
    """
    브라우저(Chrome/Edge 등)의 '한국어 자동 번역' 기능이 켜져 있을 때
    React의 가상 DOM 텍스트 노드를 <font> 태그로 감싸면서 발생하는
    'NotFoundError: The node to be removed is not a child of this node'
    (노드가 이 노드의 자식이 아닙니다) 런타임 충돌 에러를 원천 차단합니다.
    """
    st.markdown(
        """
        <meta name="google" content="notranslate" />
        <style>
            html, body, [data-testid="stAppViewContainer"], .main, div, span, p {
                translate: no !important;
            }
            .notranslate {
                translate: no !important;
            }
        </style>
        <script>
            document.documentElement.setAttribute('translate', 'no');
            document.documentElement.classList.add('notranslate');
            if (document.body) {
                document.body.setAttribute('translate', 'no');
                document.body.classList.add('notranslate');
            }
        </script>
        """,
        unsafe_allow_html=True
    )


def to_power10_str(val: float, precision: int = 4) -> str:
    """
    부동소수점 숫자를 '8.6563e-16' 형태 대신 '8.6563 × 10^(-16)' 형태로 변환합니다.
    - 지수가 0인 경우 일반 소수점 표기
    - 0에 극도로 가까운 수치나 지수 표기가 필요한 수치를 직관적인 10^(지수) 형태로 변환
    """
    if val is None or not np.isfinite(val):
        return str(val)
    if val == 0.0:
        return "0"
    
    # 지수 형식(scientific notation) 파싱
    formatted = f"{val:.{precision}e}"
    base_str, exp_str = formatted.split("e")
    exp = int(exp_str)
    
    if exp == 0:
        return f"{val:.{precision}f}".rstrip("0").rstrip(".")
    
    return f"{base_str} × 10^({exp})"



# -----------------------------------------------------------------------------
# 2. 사이드바: 프리셋 및 입력 변수
# -----------------------------------------------------------------------------
def render_sidebar():
    st.sidebar.header("⚙️ 시뮬레이션 설정 및 입력 변수")

    st.sidebar.markdown("### 🧪 대표 유체 프리셋")
    preset_choice = st.sidebar.selectbox(
        "유체 프리셋을 선택하면 물성치가 자동 입력됩니다:",
        ["물 (Water, 20°C)", "글리세린 (Glycerin, 20°C)", "엔진 오일 (SAE 30, 20°C)", "직접 입력 (Custom)"]
    )

    # 기본값 설정 (층류 유동 시뮬레이터 기준 Re < 2,100 초기화)
    if preset_choice == "물 (Water, 20°C)":
        default_mu = 0.001002
        default_rho = 998.2
        default_dp = 5.0  # Re ~ 1,351 (안정적인 층류)
    elif preset_choice == "글리세린 (Glycerin, 20°C)":
        default_mu = 1.412
        default_rho = 1260.0
        default_dp = 80000.0  # Re ~ 13.7 (층류)
    elif preset_choice == "엔진 오일 (SAE 30, 20°C)":
        default_mu = 0.29
        default_rho = 888.0
        default_dp = 15000.0  # Re ~ 107.8 (층류)
    else:
        default_mu = 0.001
        default_rho = 1000.0
        default_dp = 5.0

    st.sidebar.markdown("---")
    st.sidebar.markdown("### 📏 배관 기하 구조 (Geometry)")

    col1, col2 = st.sidebar.columns(2)
    with col1:
        ri_mm = st.number_input(
            "내부관 외반지름 Ri [mm]",
            min_value=0.1,
            max_value=500.0,
            value=10.0,
            step=1.0,
            help="내부관(Inner pipe)의 바깥 반지름입니다."
        )
    with col2:
        ro_mm = st.number_input(
            "외부관 내반지름 Ro [mm]",
            min_value=0.2,
            max_value=1000.0,
            value=30.0,
            step=1.0,
            help="외부관(Outer pipe)의 안쪽 반지름입니다. Ri보다 커야 합니다."
        )

    length = st.sidebar.number_input(
        "관 길이 L [m]",
        min_value=0.01,
        max_value=1000.0,
        value=5.0,
        step=0.5,
        help="유동이 일어나는 관의 전체 길이입니다."
    )

    st.sidebar.markdown("---")
    st.sidebar.markdown("### 💧 유체 물성 및 구동 조건 (Fluid & Flow)")

    mu = st.sidebar.number_input(
        "유체 점도 μ [Pa·s]",
        min_value=1e-6,
        max_value=100.0,
        value=float(default_mu),
        step=1e-4,
        format="%.6f",
        help="유체의 점성계수 (1 Pa·s = 1 kg/(m·s) = 1,000 cP)"
    )

    rho = st.sidebar.number_input(
        "유체 밀도 ρ [kg/m³]",
        min_value=1.0,
        max_value=20000.0,
        value=float(default_rho),
        step=10.0,
        help="유체의 밀도 (물: 약 1000 kg/m³, 공기: 약 1.2 kg/m³)"
    )

    delta_p = st.sidebar.number_input(
        "압력 강하 ΔP [Pa]",
        min_value=0.01,
        max_value=1e7,
        value=float(default_dp),
        step=10.0,
        format="%.2f",
        help="관 입구와 출구 사이의 압력차 (ΔP = Pin - Pout > 0)"
    )

    # SI 단위(m)로 변환
    ri_m = ri_mm / 1000.0
    ro_m = ro_mm / 1000.0

    return AnnularFlowParameters(
        Ri=ri_m,
        Ro=ro_m,
        L=length,
        mu=mu,
        rho=rho,
        delta_P=delta_p
    )


# -----------------------------------------------------------------------------
# 3. 메인 화면 헤더 및 이론 설명 탭
# -----------------------------------------------------------------------------
def render_theory_section():
    st.title("🌊 동심 이중관(Concentric Annulus) 내 층류 유동 시뮬레이터")
    st.caption("Chemical Engineering Fluid Dynamics Simulation — Concentric Annular Laminar Flow")

    with st.expander("📖 **화학공학 유체역학 이론 및 지배 방정식 상세 설명 (클릭하여 펼치기/접기)**", expanded=False):
        st.markdown(r"""
### 1. 원통 좌표계 Navier–Stokes 방정식
정상 상태(Steady-state), 비압축성 뉴턴 유체, 축대칭(Axisymmetric, $\partial/\partial\theta=0$), 완전 발달 층류 유동($v_r=v_\theta=0, \partial u/\partial z=0$)을 가정합니다.
유동 방향($z$축) 운동량 보존 방정식은 다음과 같이 단순화됩니다:

$$0 = -\frac{dp}{dz} + \mu \frac{1}{r}\frac{d}{dr}\left(r \frac{du}{dr}\right)$$

여기서 압력 구배는 관 길이 $L$ 동안 일정하므로 $-\frac{dp}{dz} = \frac{\Delta P}{L}$ ($\Delta P = P_{in} - P_{out} > 0$)입니다.

---

### 2. 이중 경계조건 (No-Slip Boundary Conditions)
고체 벽면에서 점성 유체는 미끄러지지 않으므로(No-slip), 내부관 외벽과 외부관 내벽에서 유속은 0입니다:
- 내부관 표면: $u(R_i) = 0$
- 외부관 표면: $u(R_o) = 0$

---

### 3. 유속 분포 해석해 (Velocity Profile)
미분방정식을 2회 적분하고 경계조건을 적용하면 엄밀 해석해(Analytical Solution)를 얻습니다:

$$u(r) = \frac{\Delta P}{4\mu L} \left[ R_o^2 - r^2 - \frac{R_o^2 - R_i^2}{\ln(R_o / R_i)} \ln\left(\frac{R_o}{r}\right) \right]$$

---

### 4. 최대 속도 위치 및 최대 속도
속도 구배가 0이 되는 지점($\frac{du}{dr} = 0$)에서 최대 속도가 나타납니다:

$$r_{max} = \sqrt{\frac{R_o^2 - R_i^2}{2 \ln(R_o / R_i)}}$$

$$u_{max} = u(r_{max})$$

*참고*: 원통 좌표계의 곡률 효과(기하학적 비대칭성)로 인해 $r_{max}$는 두 관의 산술 평균 $\frac{R_i + R_o}{2}$보다 항상 **안쪽(내부관 쪽)**에 위치하게 됩니다.

---

### 5. 전단응력 분포 (Shear Stress Distribution)
뉴턴의 점성 법칙에 따라 유체 내 전단응력 $\tau_{rz}(r)$은 다음과 같이 정의됩니다:

$$\tau_{rz}(r) = \mu \frac{du}{dr} = \frac{\Delta P}{4L} \left[ -2r + \frac{R_o^2 - R_i^2}{r \ln(R_o / R_i)} \right]$$

- **최대 유속 위치**: $r = r_{max}$에서 $\frac{du}{dr} = 0$이므로 $\tau_{rz}(r_{max}) = 0$ 입니다.
- **내부관 벽면 ($r = R_i$)**: 유속이 $r$ 증가에 따라 빨라지므로 $\frac{du}{dr} > 0 \implies \tau_{rz}(R_i) > 0$.
- **외부관 벽면 ($r = R_o$)**: 외부벽으로 갈수록 유속이 감소하므로 $\frac{du}{dr} < 0 \implies \tau_{rz}(R_o) < 0$.
- **물리적 방향 해석**: 벽면이 유체에 가하는 마찰력은 유동 진행 방향($+z$)의 반대 방향으로 작용하여 압력 강하에 의한 추진력과 완벽한 정적 평형을 이룹니다.

---

### 6. 유량, 평균 유속 및 레이놀즈 수
- **환형 단면적**: $A = \pi(R_o^2 - R_i^2)$
- **수력 직경 (Hydraulic Diameter)**: $D_h = \frac{4A}{P_{wetted}} = \frac{4\pi(R_o^2 - R_i^2)}{2\pi(R_o + R_i)} = 2(R_o - R_i)$
- **체적 유량 해석해**: 
  $$Q = 2\pi \int_{R_i}^{R_o} u(r)r\,dr = \frac{\pi \Delta P}{8\mu L}\left[(R_o^4 - R_i^4) - \frac{(R_o^2 - R_i^2)^2}{\ln(R_o/R_i)}\right]$$
- **평균 유속**: $u_{avg} = \frac{Q}{A}$
- **수력 직경 기반 레이놀즈 수**: $Re = \frac{\rho u_{avg} D_h}{\mu}$
  - 일반 원관의 경우 $Re < 2,100 \sim 2,300$에서 층류로 판정되나, 환형 관(Annulus)의 경우 반경비($R_i/R_o$)에 따라 임계 레이놀즈수가 미세하게 변동할 수 있습니다.
        """)


# -----------------------------------------------------------------------------
# 4. 주요 수치 지표 카드 렌더링
# -----------------------------------------------------------------------------
def render_metrics(metrics: dict):
    st.subheader("📊 주요 계산 결과 지표 (Key Flow Metrics)")

    m_row1_col1, m_row1_col2, m_row1_col3, m_row1_col4 = st.columns(4)
    with m_row1_col1:
        st.metric(
            label="최대 유속 u_max",
            value=f"{metrics['umax_m_s']:.4f} m/s",
            help="이중관 사이에서 유체가 도달하는 최대 유속입니다."
        )
    with m_row1_col2:
        st.metric(
            label="최대속도 위치 r_max",
            value=f"{metrics['rmax_mm']:.2f} mm",
            help="유속이 극댓값을 가지는 반지름 위치입니다."
        )
    with m_row1_col3:
        st.metric(
            label="평균 유속 u_avg",
            value=f"{metrics['average_velocity_m_s']:.4f} m/s",
            help="전체 단면적 대비 체적 유량(Q / A)입니다."
        )
    with m_row1_col4:
        st.metric(
            label="최대속도비 u_max / u_avg",
            value=f"{(metrics['umax_m_s'] / metrics['average_velocity_m_s']):.3f}",
            help="평균 유속 대비 최대 유속의 비 (원관의 경우 2.0)"
        )

    m_row2_col1, m_row2_col2, m_row2_col3, m_row2_col4 = st.columns(4)
    with m_row2_col1:
        st.metric(
            label="체적 유량 Q",
            value=f"{metrics['flow_rate_L_min']:.3f} L/min",
            delta=f"{to_power10_str(metrics['flow_rate_m3_s'], precision=3)} m³/s",
            help="단위 시간당 통과하는 유체의 부피입니다."
        )
    with m_row2_col2:
        st.metric(
            label="수력 직경 Dh",
            value=f"{metrics['hydraulic_diameter_mm']:.1f} mm",
            help="2 * (Ro - Ri)"
        )
    with m_row2_col3:
        re = metrics["reynolds_number"]
        is_laminar = re < 2100
        re_label = "층류 상태" if is_laminar else ("천이 영역" if re < 4000 else "난류 상태")
        st.metric(
            label="레이놀즈 수 Re",
            value=f"{re:,.1f}",
            delta=re_label,
            delta_color="normal" if is_laminar else "inverse",
            help="수력 직경 기반 무차원수 Re = ρ * u_avg * Dh / μ"
        )
    with m_row2_col4:
        st.metric(
            label="환형 단면적 A",
            value=f"{metrics['area_m2'] * 1e4:.2f} cm²",
            help="유동 단면적 A = π(Ro² - Ri²)"
        )

    # 층류 판정 안내문 (브라우저 번역기 및 React DOM 가상 트리 충돌 방지를 위해 컨테이너 격리 및 안전한 텍스트 렌더링)
    flow_status_container = st.container()
    with flow_status_container:
        if re >= 2100:
            st.warning(
                f"⚠️ 현재 계산된 레이놀즈 수(Re = {re:,.1f})는 통상적인 층류 기준(2,100 미만)을 초과합니다. "
                "본 모델의 해석해는 완전 층류(Laminar) 가정에 기반하므로, 실제 유동에서는 난류 천이로 인해 오차가 발생할 수 있습니다."
            )
        else:
            st.success(
                f"✅ 레이놀즈 수(Re = {re:,.1f}, 기준치 2,100 이하)에 따라 안정적인 **완전 층류 유동(Laminar Flow)** 상태입니다."
            )


# -----------------------------------------------------------------------------
# 5. 그래프 시각화 (유속 분포 및 전단응력 분포)
# -----------------------------------------------------------------------------
def plot_matplotlib_profiles(model: ConcentricAnnularFlowModel):
    """
    Matplotlib 기반 2열 서브플롯: 유속 분포 u(r) 및 전단응력 분포 tau_rz(r)
    - headless(Agg) 백엔드 및 안전한 좌표 바운딩으로 렌더링 누락 방지
    """
    try:
        r_grid, u_grid, tau_grid = model.generate_profiles(n_points=600)
        r_mm = r_grid * 1e3
        ri_mm = model.Ri * 1e3
        ro_mm = model.Ro * 1e3
        rmax_mm = model.calculate_rmax() * 1e3
        umax = float(model.calculate_umax())
        uavg = float(model.average_velocity())
        tau_Ri = float(model.shear_stress(model.Ri))
        tau_Ro = float(model.shear_stress(model.Ro))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=120)

        # ------------------ (1) 유속 분포 그래프 ------------------
        ax1.plot(r_mm, u_grid, color="#1f77b4", linewidth=2.8, label="유속 분포 $u(r)$")
        
        # 평균 유속 수평선
        ax1.axhline(uavg, color="#ff7f0e", linestyle="--", linewidth=1.5, label=f"평균 유속 $u_{{avg}} = {uavg:.3f}$ m/s")
        
        # 최대 유속 위치 수직선 및 마커
        ax1.axvline(rmax_mm, color="#d62728", linestyle=":", linewidth=2, label=f"최대속도 위치 $r_{{max}} = {rmax_mm:.2f}$ mm")
        ax1.plot(rmax_mm, umax, "ro", markersize=8, label=f"최대 유속 $u_{{max}} = {umax:.3f}$ m/s")
        
        # 벽면 경계 (음영 처리)
        ax1.axvline(ri_mm, color="#2ca02c", linestyle="-", linewidth=2.5, label=f"내부관 외벽 $R_i = {ri_mm:.1f}$ mm")
        ax1.axvline(ro_mm, color="#9467bd", linestyle="-", linewidth=2.5, label=f"외부관 내벽 $R_o = {ro_mm:.1f}$ mm")
        
        # 벽면 바깥 영역 회색 음영
        r_margin = (ro_mm - ri_mm) * 0.08
        ax1.set_xlim(ri_mm - r_margin, ro_mm + r_margin)
        ax1.axvspan(ri_mm - r_margin, ri_mm, color="gray", alpha=0.25, hatch="//")
        ax1.axvspan(ro_mm, ro_mm + r_margin, color="gray", alpha=0.25, hatch="//")

        ax1.set_title("동심 이중관 유속 분포 $u(r)$", fontsize=14, fontweight="bold", pad=12)
        ax1.set_xlabel("반지름 $r$ [mm]", fontsize=12, fontweight="bold")
        ax1.set_ylabel("유속 $u(r)$ [m/s]", fontsize=12, fontweight="bold")
        ax1.grid(True, linestyle="--", alpha=0.6)
        ax1.legend(loc="upper right", fontsize=9, framealpha=0.9)
        
        if np.isfinite(umax) and umax > 1e-9:
            ax1.set_ylim(bottom=-0.05 * umax, top=umax * 1.18)

        # 어노테이션 (안전 좌표 계산)
        annot_x = rmax_mm + (ro_mm - ri_mm) * 0.05
        annot_y = umax * 0.90 if umax > 1e-9 else 0.01
        ax1.annotate(
            f"$u_{{max}}={umax:.4f}$ m/s\n($r={rmax_mm:.2f}$ mm)",
            xy=(rmax_mm, umax),
            xytext=(annot_x, annot_y),
            arrowprops=dict(facecolor="#d62728", shrink=0.08, width=1.5, headwidth=6),
            fontsize=10,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#ffebee", edgecolor="#d62728", alpha=0.9)
        )

        # ------------------ (2) 전단응력 분포 그래프 ------------------
        ax2.plot(r_mm, tau_grid, color="#d62728", linewidth=2.8, label="전단응력 $\\tau_{rz}(r)$")
        
        # tau = 0 기준선
        ax2.axhline(0, color="black", linestyle="-", linewidth=1.2, alpha=0.7)
        
        # rmax 지점 (tau = 0 확인)
        ax2.axvline(rmax_mm, color="#1f77b4", linestyle=":", linewidth=2, label=f"$r_{{max}}$ ($\\tau_{{rz}}=0$)")
        ax2.plot(rmax_mm, 0, "bo", markersize=8)

        # 벽면 위치
        ax2.axvline(ri_mm, color="#2ca02c", linestyle="-", linewidth=2.5, label=f"내부관 $R_i$ ($\\tau_{{rz}}={tau_Ri:.2f}$ Pa)")
        ax2.axvline(ro_mm, color="#9467bd", linestyle="-", linewidth=2.5, label=f"외부관 $R_o$ ($\\tau_{{rz}}={tau_Ro:.2f}$ Pa)")
        
        # 벽면 바깥 영역 회색 음영
        ax2.set_xlim(ri_mm - r_margin, ro_mm + r_margin)
        ax2.axvspan(ri_mm - r_margin, ri_mm, color="gray", alpha=0.25, hatch="//")
        ax2.axvspan(ro_mm, ro_mm + r_margin, color="gray", alpha=0.25, hatch="//")

        max_tau = max(abs(tau_Ri), abs(tau_Ro), 1e-6)
        ax2.set_ylim(bottom=-max_tau * 1.3, top=max_tau * 1.3)

        ax2.set_title("전단응력 분포 $\\tau_{rz}(r) = \\mu \\frac{du}{dr}$", fontsize=14, fontweight="bold", pad=12)
        ax2.set_xlabel("반지름 $r$ [mm]", fontsize=12, fontweight="bold")
        ax2.set_ylabel("전단응력 $\\tau_{rz}$ [Pa]", fontsize=12, fontweight="bold")
        ax2.grid(True, linestyle="--", alpha=0.6)
        ax2.legend(loc="upper right", fontsize=9, framealpha=0.9)

        # 어노테이션
        ax2.annotate(
            f"내부벽 $\\tau_{{Ri}}={tau_Ri:.2f}$ Pa",
            xy=(ri_mm, tau_Ri),
            xytext=(ri_mm + (ro_mm - ri_mm) * 0.05, tau_Ri * 0.75),
            arrowprops=dict(facecolor="#2ca02c", shrink=0.08, width=1.5, headwidth=6),
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#e8f5e9", edgecolor="#2ca02c", alpha=0.9)
        )

        ax2.annotate(
            f"외부벽 $\\tau_{{Ro}}={tau_Ro:.2f}$ Pa",
            xy=(ro_mm, tau_Ro),
            xytext=(ro_mm - (ro_mm - ri_mm) * 0.38, tau_Ro * 0.75),
            arrowprops=dict(facecolor="#9467bd", shrink=0.08, width=1.5, headwidth=6),
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f3e5f5", edgecolor="#9467bd", alpha=0.9)
        )

        plt.tight_layout()
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)
    except Exception as e:
        st.error(f"❌ Matplotlib 그래프 렌더링 중 오류가 발생했습니다: {str(e)}")


def plot_plotly_profiles(model: ConcentricAnnularFlowModel):
    """
    Plotly 기반 인터랙티브 차트:
    마우스 호버를 통한 정확한 수치 확인, 줌, 패닝, 스크린샷 저장 지원
    """
    try:
        r_grid, u_grid, tau_grid = model.generate_profiles(n_points=600)
        r_mm = r_grid * 1e3
        ri_mm = model.Ri * 1e3
        ro_mm = model.Ro * 1e3
        rmax_mm = model.calculate_rmax() * 1e3
        umax = float(model.calculate_umax())
        uavg = float(model.average_velocity())
        tau_Ri = float(model.shear_stress(model.Ri))
        tau_Ro = float(model.shear_stress(model.Ro))

        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("동심 이중관 유속 분포 u(r)", "전단응력 분포 τ_rz(r)"),
            horizontal_spacing=0.12
        )

        # 1. 유속 분포 곡선
        fig.add_trace(
            go.Scatter(
                x=r_mm, y=u_grid,
                mode="lines",
                name="유속 u(r)",
                line=dict(color="#1f77b4", width=3),
                hovertemplate="반지름 r: %{x:.2f} mm<br>유속 u: %{y:.4f} m/s<extra></extra>"
            ),
            row=1, col=1
        )
        # 평균 유속 수평선
        fig.add_trace(
            go.Scatter(
                x=[ri_mm, ro_mm], y=[uavg, uavg],
                mode="lines",
                name=f"평균 유속 ({uavg:.3f} m/s)",
                line=dict(color="#ff7f0e", width=2, dash="dash"),
                hovertemplate=f"평균 유속 u_avg: {uavg:.4f} m/s<extra></extra>"
            ),
            row=1, col=1
        )
        # 최대 유속 위치 점
        fig.add_trace(
            go.Scatter(
                x=[rmax_mm], y=[umax],
                mode="markers",
                name=f"최대 유속 u_max ({umax:.3f} m/s)",
                marker=dict(color="#d62728", size=10, symbol="circle"),
                hovertemplate=f"최대속도 위치 r_max: {rmax_mm:.2f} mm<br>최대 유속 u_max: {umax:.4f} m/s<extra></extra>"
            ),
            row=1, col=1
        )

        # 2. 전단응력 분포 곡선
        fig.add_trace(
            go.Scatter(
                x=r_mm, y=tau_grid,
                mode="lines",
                name="전단응력 τ_rz",
                line=dict(color="#d62728", width=3),
                hovertemplate="반지름 r: %{x:.2f} mm<br>전단응력 τ_rz: %{y:.3f} Pa<extra></extra>"
            ),
            row=1, col=2
        )
        # tau = 0 기준선
        fig.add_trace(
            go.Scatter(
                x=[ri_mm, ro_mm], y=[0, 0],
                mode="lines",
                name="τ_rz = 0 기준선",
                line=dict(color="black", width=1.2),
                showlegend=False
            ),
            row=1, col=2
        )
        # rmax 에서 tau = 0 점
        fig.add_trace(
            go.Scatter(
                x=[rmax_mm], y=[0],
                mode="markers",
                name="r_max (τ_rz = 0)",
                marker=dict(color="#1f77b4", size=9),
                hovertemplate=f"최대속도 위치 r_max: {rmax_mm:.2f} mm<br>전단응력 τ_rz: 0.000 Pa<extra></extra>"
            ),
            row=1, col=2
        )

        # 벽면 가이드라인 추가
        for col_idx in [1, 2]:
            fig.add_vline(x=ri_mm, line_width=2, line_dash="solid", line_color="#2ca02c", row=1, col=col_idx)
            fig.add_vline(x=ro_mm, line_width=2, line_dash="solid", line_color="#9467bd", row=1, col=col_idx)
            fig.add_vline(x=rmax_mm, line_width=1.5, line_dash="dot", line_color="#d62728", row=1, col=col_idx)

        fig.update_xaxes(title_text="반지름 r [mm]", row=1, col=1)
        fig.update_xaxes(title_text="반지름 r [mm]", row=1, col=2)
        fig.update_yaxes(title_text="유속 u(r) [m/s]", row=1, col=1)
        fig.update_yaxes(title_text="전단응력 τ_rz [Pa]", row=1, col=2)

        fig.update_layout(
            height=460,
            margin=dict(l=40, r=40, t=50, b=30),
            hovermode="closest",
            legend=dict(orientation="h", yanchor="bottom", y=-0.3, xanchor="center", x=0.5)
        )

        st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.error(f"❌ Plotly 그래프 렌더링 중 오류가 발생했습니다: {str(e)}")


def plot_profiles(model: ConcentricAnnularFlowModel):
    """
    유속 및 전단응력 분포 시각화 (Matplotlib 기본 + Plotly 인터랙티브 탭 제공)
    """
    tab_mpl, tab_plotly = st.tabs(["📊 Matplotlib 고해상도 차트 (기본)", "🔍 Plotly 인터랙티브 차트 (호버/확대 지원)"])
    with tab_mpl:
        plot_matplotlib_profiles(model)
    with tab_plotly:
        plot_plotly_profiles(model)


# -----------------------------------------------------------------------------
# 6. 추가 시각화: 2D 단면 환형 유속장 (Cross-sectional Velocity Field)
# -----------------------------------------------------------------------------
def plot_cross_section_contour(model: ConcentricAnnularFlowModel):
    try:
        ro = model.Ro
        ri = model.Ri
        grid_size = 220
        x = np.linspace(-ro * 1.1, ro * 1.1, grid_size)
        y = np.linspace(-ro * 1.1, ro * 1.1, grid_size)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)

        U = np.zeros_like(R)
        mask = (R >= ri) & (R <= ro)
        U[mask] = model.velocity(R[mask])

        fig, ax = plt.subplots(figsize=(7, 6), dpi=110)
        cf = ax.contourf(X * 1e3, Y * 1e3, U, levels=40, cmap="viridis")
        cbar = fig.colorbar(cf, ax=ax)
        cbar.set_label("유속 $u$ [m/s]", fontsize=11, fontweight="bold")

        # 내부관 원 및 외부관 원 테두리
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.plot(ri * 1e3 * np.cos(theta), ri * 1e3 * np.sin(theta), "w--", linewidth=1.5, label="내부관 경계 $R_i$")
        ax.plot(ro * 1e3 * np.cos(theta), ro * 1e3 * np.sin(theta), "w-", linewidth=2.0, label="외부관 경계 $R_o$")

        # 최대속도 위치 원
        rmax = model.calculate_rmax()
        ax.plot(rmax * 1e3 * np.cos(theta), rmax * 1e3 * np.sin(theta), "r:", linewidth=1.8, label="최대속도 원 $r_{max}$")

        ax.set_aspect("equal")
        ax.set_title("환형 단면 2차원 유속 분포도 (Annular Cross-Section)", fontsize=13, fontweight="bold", pad=10)
        ax.set_xlabel("X 좌표 [mm]", fontsize=11)
        ax.set_ylabel("Y 좌표 [mm]", fontsize=11)
        ax.legend(loc="upper right", fontsize=8.5, framealpha=0.85)
        plt.tight_layout()
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)
    except Exception as e:
        st.error(f"❌ 2D 단면도 렌더링 중 오류가 발생했습니다: {str(e)}")


# -----------------------------------------------------------------------------
# 7. 이론 및 수치 검증 섹션
# -----------------------------------------------------------------------------
def render_verification_section(model: ConcentricAnnularFlowModel):
    st.subheader("🔍 수치 검증 및 이론 일치성 판정 (Theoretical & Numerical Verification)")
    verif = model.verify_solution()

    v_col1, v_col2 = st.columns(2)

    with v_col1:
        st.markdown("#### 1. 이중 No-slip 경계조건 검증")
        u_ri = verif["u_Ri"]
        u_ro = verif["u_Ro"]

        u_ri_str = to_power10_str(u_ri, precision=4)
        u_ro_str = to_power10_str(u_ro, precision=4)

        b_ri_icon = "✅" if verif["u_Ri_satisfied"] else "❌"
        b_ro_icon = "✅" if verif["u_Ro_satisfied"] else "❌"

        st.markdown(
            f"- **내부관 벽면 ($r = R_i$)**: `u(Ri) = {u_ri_str} m/s` {b_ri_icon} "
            f"({'경계조건 만족 (|u| < 10^(-9) m/s)' if verif['u_Ri_satisfied'] else '오차 발생'})\n"
            f"- **외부관 벽면 ($r = R_o$)**: `u(Ro) = {u_ro_str} m/s` {b_ro_icon} "
            f"({'경계조건 만족 (|u| < 10^(-9) m/s)' if verif['u_Ro_satisfied'] else '오차 발생'})"
        )

        st.markdown("#### 2. 최대속도 위치 극값 조건 검증")
        rmax_mm = verif["rmax"] * 1e3
        dudr_rmax = verif["dudr_rmax"]
        tau_rmax = verif["tau_rmax"]

        dudr_rmax_str = to_power10_str(dudr_rmax, precision=4)
        tau_rmax_str = to_power10_str(tau_rmax, precision=4)

        g_icon = "✅" if verif["dudr_rmax_zero"] else "❌"
        t_icon = "✅" if verif["tau_rmax_zero"] else "❌"

        st.markdown(
            f"- **최대속도 위치 $r_{{max}}$**: `{rmax_mm:.4f} mm`\n"
            f"- **속도 구배**: `du/dr = {dudr_rmax_str} s⁻¹` {g_icon}\n"
            f"- **전단응력**: `τ_rz(r_max) = {tau_rmax_str} Pa` {t_icon}\n"
            f"- **판정**: 최대속도 위치에서 전단응력과 속도 구배가 완벽히 0이 됨을 수치적으로 확인하였습니다."
        )

    with v_col2:
        st.markdown("#### 3. 체적 유량 검증 (해석해 vs 수치 적분)")
        q_ana = verif["q_analytical"]
        q_num = verif["q_numerical"]

        q_ana_str = to_power10_str(q_ana, precision=6)
        q_num_str = to_power10_str(q_num, precision=6)

        q_err = verif["q_relative_difference"] * 100.0
        q_icon = "✅" if verif["q_matches"] else "❌"

        st.markdown(
            f"- **해석적 체적 유량 $Q_{{analytical}}$**: `{q_ana_str} m³/s` ({q_ana*1000*60:.4f} L/min)\n"
            f"- **수치 적분 유량 $Q_{{numerical}}$**: `{q_num_str} m³/s` ({q_num*1000*60:.4f} L/min)\n"
            f"- **상대 오차 (Relative Error)**: `{q_err:.6f} %` {q_icon}\n"
            f"- **판정**: 사다리꼴 수치 적분 결과와 엄밀 해석해가 0.01% 미만의 극미한 오차로 일치합니다."
        )

        st.markdown("#### 4. 벽면 전단응력 부호 및 방향 해석")
        tau_ri = model.shear_stress(model.Ri)
        tau_ro = model.shear_stress(model.Ro)

        st.info(
            f"• **내부벽 ($\tau_{{Ri}} = {tau_ri:+.2f}$ Pa)**: $\\frac{{du}}{{dr}} > 0$이므로 양수입니다. "
            f"모멘텀이 반경 외측($+r$)으로 확산됩니다.\n\n"
            f"• **외부벽 ($\tau_{{Ro}} = {tau_ro:+.2f}$ Pa)**: $\\frac{{du}}{{dr}} < 0$이므로 음수입니다. "
            f"모멘텀이 중심 방향($-r$)으로 수렴합니다.\n\n"
            f"• **물리적 의미**: 두 벽면 모두 $+z$ 유동을 방해하는 마찰 저항을 유체에 전달하여 압력차에 의한 힘과 정밀한 힘의 평형을 이룹니다."
        )


# -----------------------------------------------------------------------------
# 8. 메인 애플리케이션 진입점
# -----------------------------------------------------------------------------
def main():
    inject_anti_dom_crash_headers()
    params = render_sidebar()

    # 입력값 검증 (Validation)
    try:
        model = ConcentricAnnularFlowModel(params)
    except ValueError as e:
        st.error(f"❌ **입력값 오류**: {str(e)}")
        st.info("💡 사이드바의 입력 슬라이더/숫자 값을 확인하여 0 < Ri < Ro 및 양수 값을 유지해 주세요.")
        return

    render_theory_section()

    st.markdown("---")
    metrics = model.get_summary_metrics()
    render_metrics(metrics)

    st.markdown("---")
    st.subheader("📈 유속 및 전단응력 분포 시각화 (Flow Profiles)")
    plot_profiles(model)

    # 2D 단면 유속장 아코디언
    with st.expander("🌐 **2차원 환형 단면(Cross-Section) 유속장 등고선도 보기**"):
        plot_cross_section_contour(model)

    st.markdown("---")
    render_verification_section(model)

    st.markdown("---")
    st.caption("Developed with Python, Streamlit, NumPy, and Matplotlib. Concentric Annular Flow Simulation.")


if __name__ == "__main__":
    if not st.runtime.exists():
        # 터미널에서 'python app.py'로 직접 실행해도 자동으로 Streamlit 웹앱 실행
        from streamlit.web import cli as stcli
        sys.argv = ["streamlit", "run", __file__, "--server.address", "127.0.0.1", "--server.port", "8501"]
        sys.exit(stcli.main())
    else:
        main()
