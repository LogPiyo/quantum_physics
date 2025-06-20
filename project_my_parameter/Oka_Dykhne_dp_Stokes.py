# %%[markdown]
# これは解析的に求めた占有確率において，Stokes位相をフィッティングパラメータとしたときの
# Stokes位相のエネルギー依存性を可視化するプログラムです。

# %%
# multiple-passage TLZ modelにおける断熱状態の占有確率を数値計算結果と理論値で比較する
import _pathmagic # noqa
import math
import cmath
import scipy
import numpy as np
import matplotlib.pyplot as plt

from my_module.function import q, adia_eng, adia_param, to_LZ
from my_module.calculator import calculate_occupation_probability
from scipy.integrate import quad


# parameter
eps_0 = -50
D_y_val = np.linspace(-100, 100, 50)  # twist strength
# D_y_val = D_y_val[(D_y_val >= 5) | (D_y_val <= -5)]
D_z = 4  # minimal energy gap
F = -1  # sweep speed (should not change)(default value: -1)(時間反転させないため)
t_i = -math.pi / abs(F)  # initial time
t_f = math.pi / abs(F)  # final time

tp_1 = -math.pi / (2 * abs(F))  # first transition time
tp_2 = math.pi / (2 * abs(F))  # second transition time

# constant
h = 1  # Dirac constant (should not change, initial value: 1)
n = 500  # step
OP_list = []  # ocupation probability
Stokes_val = []
Stokes_val_thr_TLZ = []
Stokes_val_thr_LZ = []
t_eval = np.linspace(t_i, t_f, n)  # time


def Hc(t, component):
    """
    Hamiltonianの設定(複素数対応)

    Args:
        t (float): time
        component (string): 成分

    Returns:
        float: 時刻tにおけるcomponentで指定した成分を返す。
    """
    H = {}

    H['x'] = -eps_0 * cmath.cos(q(t, F))
    H['y'] = -0.125 * (4 * D_y / eps_0**2) * eps_0**2 * cmath.sin(2 * q(t, F))**2
    H['z'] = D_z * cmath.sin(q(t, F))
    H['x_dot'] = eps_0 * cmath.sin(q(t, F))
    H['y_dot'] = -0.125 * (4 * D_y / eps_0**2) * eps_0**2 * 4 * cmath.sin(2 * q(t, F)) * cmath.cos(2 * q(t, F))
    H['z_dot'] = D_z * cmath.cos(q(t, F))

    return H[component]


def Re_E(t):
    """
    define real part of adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: adiabatic energy
    """
    Integrand = adia_eng(tp_1 + 1j * t, to_LZ(Hc, F))
    return Integrand.real


def Im_E_1(t):
    Integrand = adia_eng(tp_1 + 1j * t, to_LZ(Hc, F))
    return Integrand.imag


def Im_E_2(t):
    Integrand = adia_eng(tp_2 + 1j * t, to_LZ(Hc, F))
    return Integrand.imag


def E_3(t):
    Integrand = adia_eng(t, to_LZ(Hc, F))
    return Integrand.real


def Stokes_phase(D_y):
    delta = adia_param(eps_0, F, D_z, (4 * D_y / eps_0**2))

    term1 = math.pi / 4
    term2 = delta * (math.log(delta) - 1)
    term3 = cmath.phase(scipy.special.gamma(1 - 1j * delta))
    return term1 + term2 + term3


for initial_D_y in [-0.1, 0.1]:
    # 基準
    prev_ans_phi_s = Stokes_phase(initial_D_y)

    if initial_D_y < 0:  # start from v = -5
        # v_val_find_Stokes is [-2, -4, -6, ...]
        D_y_val_find_Stokes = D_y_val[D_y_val <= 0]
        D_y_val_find_Stokes = D_y_val_find_Stokes[:: -1]

    else:  # start from v = 5
        # v_val_find_Stokes is [2, 4, 6, ...]
        D_y_val_find_Stokes = D_y_val[D_y_val > 0]

    for D_y in D_y_val_find_Stokes:
        ans_phi_s = 0
        min_error = 10
        for d_phi_s in np.linspace(-0.1, 0.1, 10):
            phi_s = prev_ans_phi_s + d_phi_s  # Stokes phase(弧度法)
            TLZ = -math.pi * (D_z - (4 * D_y / eps_0**2)*eps_0*F/4)**2 / (abs(eps_0) * abs(F))
            # １回目の遷移がOkaモデルと全体の符号が反転している場合は分子の第２項の符号をマイナスにする
            zero_approx = abs(D_z - (4 * D_y / eps_0**2)*(eps_0)*F/4) / (abs(eps_0) * (-F))

            # integral of Re_E
            ll_Re_E = 0  # lower limit
            ul_Re_E = zero_approx  # upper limit
            TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)
            TP *= -4 * (-F) / abs(F)

            # integral of Im_E_1
            ll_Im_E_1 = 0  # lower limit
            ul_Im_E_1 = zero_approx  # upper limit
            phase_term1, _ = quad(Im_E_1, ll_Im_E_1, ul_Im_E_1)
            phase_term1 *= (-F) / abs(F)

            # integral of Im_E_2
            ll_Im_E_2 = 0  # lower limit
            ul_Im_E_2 = zero_approx  # upper limit
            phase_term2, _ = quad(Im_E_2, ll_Im_E_2, ul_Im_E_2)
            phase_term2 *= (-F) / abs(F)

            # integral of Im_E_3
            ll_E_3 = tp_1  # lower limit
            ul_E_3 = tp_2  # upper limit
            phase_term3, _ = quad(E_3, ll_E_3, ul_E_3)
            phase_term3 *= (-F) / abs(F)

            OP_array = calculate_occupation_probability(Hc, t_i, t_f, n)
            OP_list += OP_array.tolist()

            # 終時間における状態0の占有確率
            phase = phase_term2 - phase_term1 + phase_term3
            P_f_HS = (4 * math.exp(TP) * (1 - math.exp(TP))
                      * math.cos(phi_s + phase)**2)
            P_TP = math.exp(TP)
            P_num = OP_list[-1]

            err = abs(P_f_HS - P_num)
            if err < min_error:
                min_error = err
                ans_phi_s = phi_s

        print("D_y = ", D_y, "Stokes = ", ans_phi_s, "error = ", min_error)
        prev_ans_phi_s = ans_phi_s  # 1つ前の点を基準にする
        Stokes_val.append(ans_phi_s)

    if initial_D_y < 0:
        Stokes_val = Stokes_val[:: -1]
print("completed")

for D_y in D_y_val:
    Stokes_val_thr_TLZ.append(Stokes_phase(D_y))

for D_y in D_y_val:
    D_y = 0
    Stokes_val_thr_LZ.append(Stokes_phase(D_y))

# %%
# グラフの描画
plt.plot(D_y_val, Stokes_val, linestyle="None", marker="x", label=r"numerical")
plt.plot(D_y_val, Stokes_val_thr_TLZ, label=r"theoretical (TLZ)", color="tab:green")
plt.plot(D_y_val, Stokes_val_thr_LZ, label=r"theoretical (LZ)", color="tab:orange")
plt.xlabel(r"$\Delta_y$")
plt.ylabel(r"Stokes phase $\varphi_s$")
plt.title(rf"$\epsilon_0 = {eps_0}, \Delta_z = {D_z}, \omega = {-F}$")
plt.legend()
plt.ylim(-0.1, 1.1)
plt.show()
