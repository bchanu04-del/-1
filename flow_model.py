"""
flow_model.py
==============
동심 이중관(Concentric Annulus) 내 층류 유동 해석 및 시뮬레이션 모델 모듈.

화학공학 유체역학 이론에 기반하여 다음 항목을 계산합니다:
1. 입력 변수 유효성 검사 (Validation)
2. 해석적 유속 분포 u(r) 및 전단응력 분포 tau_rz(r)
3. 최대 유속 위치 r_max 및 최대 유속 u_max
4. 단면적, 수력 직경, 체적 유량(해석해 및 수치 적분), 평균 유속
5. 레이놀즈 수(Reynolds Number, Re) 및 층류 판정 가이드
6. 경계조건 및 극값(du/dr = 0) 수치 검증
"""

from dataclasses import dataclass
from typing import Dict, Tuple, Any
import numpy as np


@dataclass
class AnnularFlowParameters:
    """
    동심 이중관 유동 해석에 필요한 입력 매개변수 데이터 클래스 (SI 단위).

    Attributes:
        Ri (float): 내부관 반지름 (m), 0 < Ri < Ro
        Ro (float): 외부관 내벽 반지름 (m), Ro > Ri
        L (float): 관 길이 (m), L > 0
        mu (float): 유체 점도 (Pa·s 또는 kg/(m·s)), mu > 0
        rho (float): 유체 밀도 (kg/m^3), rho > 0
        delta_P (float): 압력 강하 (Pa, Pin - Pout), delta_P > 0
    """
    Ri: float
    Ro: float
    L: float
    mu: float
    rho: float
    delta_P: float

    def validate(self) -> None:
        """
        입력 매개변수의 물리적 타당성을 검사합니다.
        조건을 만족하지 않는 경우 명확한 ValueError를 발생시킵니다.
        """
        if self.Ri <= 0:
            raise ValueError(f"내부관 반지름(Ri)은 0보다 커야 합니다. (입력값: {self.Ri} m)")
        if self.Ro <= 0:
            raise ValueError(f"외부관 내벽 반지름(Ro)은 0보다 커야 합니다. (입력값: {self.Ro} m)")
        if self.Ri >= self.Ro:
            raise ValueError(
                f"내부관 반지름(Ri={self.Ri} m)은 외부관 내벽 반지름(Ro={self.Ro} m)보다 작아야 합니다. (0 < Ri < Ro)"
            )
        if self.L <= 0:
            raise ValueError(f"관 길이(L)는 0보다 커야 합니다. (입력값: {self.L} m)")
        if self.mu <= 0:
            raise ValueError(f"유체 점도(μ)는 0보다 커야 합니다. (입력값: {self.mu} Pa·s)")
        if self.rho <= 0:
            raise ValueError(f"유체 밀도(ρ)는 0보다 커야 합니다. (입력값: {self.rho} kg/m³)")
        if self.delta_P <= 0:
            raise ValueError(f"압력 강하(ΔP)는 0보다 커야 합니다. (입력값: {self.delta_P} Pa)")


class ConcentricAnnularFlowModel:
    """
    동심 이중관 층류 유동 해석 모델 클래스.
    
    정상 상태(steady-state), 비압축성 뉴턴 유체, 축대칭(axisymmetric),
    완전 발달(fully developed) 층류 유동을 가정합니다.
    """

    def __init__(self, params: AnnularFlowParameters):
        """
        모델 초기화 및 매개변수 유효성 검사
        """
        params.validate()
        self.params = params

        # 계산 편의를 위한 기본 파라미터 캐싱
        self.Ri = params.Ri
        self.Ro = params.Ro
        self.L = params.L
        self.mu = params.mu
        self.rho = params.rho
        self.delta_P = params.delta_P

        # 자주 사용되는 기하학적 상수
        self.ln_ratio = np.log(self.Ro / self.Ri)
        self.Ro2_minus_Ri2 = self.Ro**2 - self.Ri**2
        self.factor = self.delta_P / (4.0 * self.mu * self.L)

    def velocity(self, r: np.ndarray | float) -> np.ndarray | float:
        """
        임의의 반지름 위치 r에서 유속 u(r) 계산 (m/s).

        이론식:
            u(r) = (ΔP / 4μL) * [ Ro^2 - r^2 - ((Ro^2 - Ri^2) / ln(Ro/Ri)) * ln(Ro/r) ]
        """
        r_arr = np.asarray(r, dtype=float)
        # Ri <= r <= Ro 영역 밖의 값에 대해 안전하게 계산
        term1 = self.Ro**2 - r_arr**2
        term2 = (self.Ro2_minus_Ri2 / self.ln_ratio) * np.log(self.Ro / r_arr)
        u = self.factor * (term1 - term2)
        return float(u) if np.ndim(r) == 0 else u

    def velocity_gradient(self, r: np.ndarray | float) -> np.ndarray | float:
        """
        반지름 위치 r에서의 속도 구배 du/dr 계산 (1/s).

        이론식:
            du/dr = (ΔP / 4μL) * [ -2r + (Ro^2 - Ri^2) / (r * ln(Ro/Ri)) ]
        """
        r_arr = np.asarray(r, dtype=float)
        grad = self.factor * (-2.0 * r_arr + self.Ro2_minus_Ri2 / (r_arr * self.ln_ratio))
        return float(grad) if np.ndim(r) == 0 else grad

    def shear_stress(self, r: np.ndarray | float) -> np.ndarray | float:
        """
        반지름 위치 r에서의 전단응력 tau_rz(r) 계산 (Pa).

        이론식:
            tau_rz(r) = μ * (du/dr)
                      = (ΔP / 4L) * [ -2r + (Ro^2 - Ri^2) / (r * ln(Ro/Ri)) ]
        """
        return self.mu * self.velocity_gradient(r)

    def calculate_rmax(self) -> float:
        """
        최대 유속 위치 r_max 계산 (m).
        속도 구배 du/dr = 0 조건에서 유도됨:
            r_max = sqrt( (Ro^2 - Ri^2) / (2 * ln(Ro/Ri)) )
        """
        return float(np.sqrt(self.Ro2_minus_Ri2 / (2.0 * self.ln_ratio)))

    def calculate_umax(self) -> float:
        """
        최대 유속 u_max 계산 (m/s).
        """
        rmax = self.calculate_rmax()
        return float(self.velocity(rmax))

    def cross_sectional_area(self) -> float:
        """
        환형 유동 단면적 A = π(Ro^2 - Ri^2) (m^2)
        """
        return float(np.pi * self.Ro2_minus_Ri2)

    def hydraulic_diameter(self) -> float:
        """
        수력 직경 Dh = 4 * A / P_wetted = 2 * (Ro - Ri) (m)
        """
        return float(2.0 * (self.Ro - self.Ri))

    def flow_rate_analytical(self) -> float:
        """
        해석적 체적 유량 Q_analytical 계산 (m^3/s).

        이론식:
            Q = (π ΔP / 8μL) * [ Ro^4 - Ri^4 - (Ro^2 - Ri^2)^2 / ln(Ro/Ri) ]
        """
        term1 = self.Ro**4 - self.Ri**4
        term2 = (self.Ro2_minus_Ri2**2) / self.ln_ratio
        return float((np.pi * self.delta_P / (8.0 * self.mu * self.L)) * (term1 - term2))

    def flow_rate_numerical(self, n_points: int = 2000) -> float:
        """
        사다리꼴 수치 적분을 통한 체적 유량 Q_numerical 계산 (m^3/s).
            Q = 2π ∫[Ri, Ro] u(r) * r dr
        """
        r_grid = np.linspace(self.Ri, self.Ro, n_points)
        u_grid = self.velocity(r_grid)
        integrand = u_grid * r_grid
        integral = np.trapezoid(integrand, r_grid) if hasattr(np, 'trapezoid') else np.trapz(integrand, r_grid)
        return float(2.0 * np.pi * integral)

    def average_velocity(self) -> float:
        """
        평균 유속 u_avg = Q / A (m/s)
        """
        return float(self.flow_rate_analytical() / self.cross_sectional_area())

    def reynolds_number(self) -> float:
        """
        수력 직경 기반 레이놀즈 수 Re = rho * u_avg * Dh / mu
        """
        u_avg = self.average_velocity()
        Dh = self.hydraulic_diameter()
        return float(self.rho * u_avg * Dh / self.mu)

    def get_summary_metrics(self) -> Dict[str, Any]:
        """
        전체 주요 계산 지표 요약 딕셔너리 반환
        """
        rmax = self.calculate_rmax()
        umax = self.calculate_umax()
        area = self.cross_sectional_area()
        dh = self.hydraulic_diameter()
        q_ana = self.flow_rate_analytical()
        q_num = self.flow_rate_numerical()
        u_avg = self.average_velocity()
        re = self.reynolds_number()

        # 벽면 전단응력
        tau_Ri = float(self.shear_stress(self.Ri))
        tau_Ro = float(self.shear_stress(self.Ro))

        return {
            "rmax_m": rmax,
            "rmax_mm": rmax * 1e3,
            "umax_m_s": umax,
            "area_m2": area,
            "hydraulic_diameter_m": dh,
            "hydraulic_diameter_mm": dh * 1e3,
            "flow_rate_m3_s": q_ana,
            "flow_rate_L_s": q_ana * 1e3,
            "flow_rate_L_min": q_ana * 1e3 * 60.0,
            "flow_rate_numerical_m3_s": q_num,
            "average_velocity_m_s": u_avg,
            "reynolds_number": re,
            "tau_inner_Pa": tau_Ri,
            "tau_outer_Pa": tau_Ro,
        }

    def verify_solution(self) -> Dict[str, Any]:
        """
        이중 경계조건, 최대속도 위치의 속도 구배 및 전단응력,
        해석해와 수치적분 유량의 일치 여부를 검증합니다.
        """
        u_Ri = float(self.velocity(self.Ri))
        u_Ro = float(self.velocity(self.Ro))
        rmax = self.calculate_rmax()
        dudr_rmax = float(self.velocity_gradient(rmax))
        tau_rmax = float(self.shear_stress(rmax))

        q_ana = self.flow_rate_analytical()
        q_num = self.flow_rate_numerical()
        q_rel_diff = abs(q_ana - q_num) / q_ana if q_ana != 0 else 0.0

        tol = 1e-9

        return {
            "u_Ri": u_Ri,
            "u_Ri_satisfied": abs(u_Ri) < tol,
            "u_Ro": u_Ro,
            "u_Ro_satisfied": abs(u_Ro) < tol,
            "rmax": rmax,
            "dudr_rmax": dudr_rmax,
            "dudr_rmax_zero": abs(dudr_rmax) < tol,
            "tau_rmax": tau_rmax,
            "tau_rmax_zero": abs(tau_rmax) < tol,
            "q_analytical": q_ana,
            "q_numerical": q_num,
            "q_relative_difference": q_rel_diff,
            "q_matches": q_rel_diff < 1e-4,
        }

    def generate_profiles(self, n_points: int = 500) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        그래프 시각화를 위한 반지름 격자, 유속, 전단응력 배열 생성
        (반지름 단위: mm, 유속: m/s, 전단응력: Pa)
        """
        r_grid = np.linspace(self.Ri, self.Ro, n_points)
        u_grid = self.velocity(r_grid)
        tau_grid = self.shear_stress(r_grid)
        return r_grid, u_grid, tau_grid
