"""
test_flow_model.py
==================
동심 이중관 층류 유동 해석 모델(flow_model.py)에 대한 pytest 단위 테스트 모음.
"""

import pytest
import numpy as np
from flow_model import AnnularFlowParameters, ConcentricAnnularFlowModel


@pytest.fixture
def standard_params():
    """표준 테스트용 유동 조건 (물 기준)"""
    return AnnularFlowParameters(
        Ri=0.01,       # 10 mm
        Ro=0.03,       # 30 mm
        L=5.0,         # 5 m
        mu=0.001002,   # 20°C 물 점도 ~ 1.002 mPa·s
        rho=998.2,     # 20°C 물 밀도 ~ 998.2 kg/m³
        delta_P=120.0  # 120 Pa
    )


@pytest.fixture
def standard_model(standard_params):
    return ConcentricAnnularFlowModel(standard_params)


class TestAnnularFlowParameters:
    def test_valid_parameters(self, standard_params):
        """정상 매개변수 유효성 검사 통과 여부"""
        standard_params.validate()

    def test_invalid_radii_ri_greater_than_ro(self):
        """Ri >= Ro 에러 검증"""
        params = AnnularFlowParameters(Ri=0.04, Ro=0.02, L=1.0, mu=0.001, rho=1000.0, delta_P=100.0)
        with pytest.raises(ValueError, match="작아야 합니다"):
            params.validate()

    def test_invalid_radii_equal(self):
        """Ri == Ro 에러 검증"""
        params = AnnularFlowParameters(Ri=0.02, Ro=0.02, L=1.0, mu=0.001, rho=1000.0, delta_P=100.0)
        with pytest.raises(ValueError):
            params.validate()

    def test_invalid_negative_or_zero_values(self):
        """0 이하 값(Ri, Ro, L, mu, rho, delta_P) 에러 검증"""
        # Ri <= 0
        with pytest.raises(ValueError):
            AnnularFlowParameters(Ri=0.0, Ro=0.02, L=1.0, mu=0.001, rho=1000.0, delta_P=100.0).validate()
        
        # Ro <= 0
        with pytest.raises(ValueError):
            AnnularFlowParameters(Ri=-0.01, Ro=-0.02, L=1.0, mu=0.001, rho=1000.0, delta_P=100.0).validate()

        # L <= 0
        with pytest.raises(ValueError):
            AnnularFlowParameters(Ri=0.01, Ro=0.02, L=0.0, mu=0.001, rho=1000.0, delta_P=100.0).validate()

        # mu <= 0
        with pytest.raises(ValueError):
            AnnularFlowParameters(Ri=0.01, Ro=0.02, L=1.0, mu=-0.001, rho=1000.0, delta_P=100.0).validate()

        # rho <= 0
        with pytest.raises(ValueError):
            AnnularFlowParameters(Ri=0.01, Ro=0.02, L=1.0, mu=0.001, rho=0.0, delta_P=100.0).validate()

        # delta_P <= 0
        with pytest.raises(ValueError):
            AnnularFlowParameters(Ri=0.01, Ro=0.02, L=1.0, mu=0.001, rho=1000.0, delta_P=-50.0).validate()


class TestBoundaryConditions:
    def test_no_slip_boundary_conditions(self, standard_model):
        """내부관(Ri)과 외부관(Ro) 벽면에서의 유속 u = 0 (No-slip 경계조건) 검증"""
        u_Ri = standard_model.velocity(standard_model.Ri)
        u_Ro = standard_model.velocity(standard_model.Ro)

        assert abs(u_Ri) < 1e-12, f"u(Ri)가 0이 아닙니다: {u_Ri}"
        assert abs(u_Ro) < 1e-12, f"u(Ro)가 0이 아닙니다: {u_Ro}"


class TestMaximumVelocity:
    def test_rmax_location_and_velocity_gradient(self, standard_model):
        """최대 유속 위치 r_max가 관 사이에 위치하며 du/dr = 0, tau_rz = 0 인지 검증"""
        rmax = standard_model.calculate_rmax()
        assert standard_model.Ri < rmax < standard_model.Ro, "rmax는 Ri와 Ro 사이에 있어야 합니다."

        # du/dr(rmax) == 0
        dudr_rmax = standard_model.velocity_gradient(rmax)
        assert abs(dudr_rmax) < 1e-12, f"du/dr(rmax) != 0: {dudr_rmax}"

        # tau_rz(rmax) == 0
        tau_rmax = standard_model.shear_stress(rmax)
        assert abs(tau_rmax) < 1e-12, f"tau_rz(rmax) != 0: {tau_rmax}"

    def test_umax_is_global_maximum(self, standard_model):
        """u(rmax)가 임의의 모든 r에서의 유속보다 큰지 검증"""
        rmax = standard_model.calculate_rmax()
        umax = standard_model.calculate_umax()

        r_test = np.linspace(standard_model.Ri, standard_model.Ro, 100)
        u_test = standard_model.velocity(r_test)

        assert np.all(u_test <= umax + 1e-12), "u(rmax)가 최대 유속이 아닙니다."
        assert umax > 0, "유속은 양수여야 합니다."


class TestShearStress:
    def test_wall_shear_stress_signs(self, standard_model):
        """
        벽면 전단응력의 부호 및 방향 검증:
        - 내부관 벽면(r=Ri): du/dr > 0 이므로 tau_rz > 0
        - 외부관 벽면(r=Ro): du/dr < 0 이므로 tau_rz < 0
        """
        tau_Ri = standard_model.shear_stress(standard_model.Ri)
        tau_Ro = standard_model.shear_stress(standard_model.Ro)

        assert tau_Ri > 0, f"내부관 벽면 전단응력은 양수여야 합니다 (실제: {tau_Ri})"
        assert tau_Ro < 0, f"외부관 벽면 전단응력은 음수여야 합니다 (실제: {tau_Ro})"


class TestFlowRateAndReynoldsNumber:
    def test_flow_rate_analytical_vs_numerical(self, standard_model):
        """해석적 유량과 수치 적분(사다리꼴 공식) 유량의 일치성 검증 (오차 0.05% 이내)"""
        q_ana = standard_model.flow_rate_analytical()
        q_num = standard_model.flow_rate_numerical(n_points=5000)

        rel_error = abs(q_ana - q_num) / q_ana
        assert rel_error < 0.0005, f"유량 해석해와 수치해 오차가 큽니다: {rel_error*100:.4f}%"

    def test_hydraulic_diameter_and_area(self, standard_model):
        """단면적 및 수력 직경 계산 검증"""
        area = standard_model.cross_sectional_area()
        dh = standard_model.hydraulic_diameter()

        expected_area = np.pi * (standard_model.Ro**2 - standard_model.Ri**2)
        expected_dh = 2.0 * (standard_model.Ro - standard_model.Ri)

        assert np.isclose(area, expected_area)
        assert np.isclose(dh, expected_dh)

    def test_reynolds_number_calculation(self, standard_model):
        """레이놀즈 수 계산 및 타당성 검증"""
        re = standard_model.reynolds_number()
        assert re > 0, "레이놀즈 수는 양수여야 합니다."
        u_avg = standard_model.average_velocity()
        dh = standard_model.hydraulic_diameter()
        expected_re = standard_model.rho * u_avg * dh / standard_model.mu
        assert np.isclose(re, expected_re)


class TestVerificationHelper:
    def test_verify_solution_structure(self, standard_model):
        """verify_solution() 함수가 모든 검증을 정상적으로 통과하는지 확인"""
        verif = standard_model.verify_solution()
        assert verif["u_Ri_satisfied"] is True
        assert verif["u_Ro_satisfied"] is True
        assert verif["dudr_rmax_zero"] is True
        assert verif["tau_rmax_zero"] is True
        assert verif["q_matches"] is True
