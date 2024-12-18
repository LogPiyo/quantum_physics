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

from my_module.function import q, adia_eng, func_psi_module, adia_param, eig_vec
from scipy.integrate import solve_ivp, quad


# parameter
eps_0_val = np.linspace(-100, 100, 50)  # energy slope
eps_0_val = eps_0_val[(eps_0_val >= 5) | (eps_0_val <= -5)]
D_z = 4  # minimal energy gap
k = 0.1  # geodesic curvature
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


def Hc(t, component, real=False):
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
    H['y'] = -0.125 * k * eps_0**2 * cmath.sin(2 * q(t, F))**2
    H['z'] = D_z * cmath.sin(q(t, F))
    H['x_dot'] = eps_0 * cmath.sin(q(t, F))
    H['y_dot'] = -0.125 * k * eps_0**2 * 4 * cmath.sin(2 * q(t, F)) * cmath.cos(2 * q(t, F))
    H['z_dot'] = D_z * cmath.cos(q(t, F))

    if real:
        return H[component].real
    else:
        return H[component]


def Re_E(t):
    """
    define real part of adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: adiabatic energy
    """
    Integrand = adia_eng(tp_1 + 1j * t, Hc, ut=True, F=F)
    return Integrand.real


def Im_E_1(t):
    Integrand = adia_eng(tp_1 + 1j * t, Hc, ut=True, F=F)
    return Integrand.imag


def Im_E_2(t):
    Integrand = adia_eng(tp_2 + 1j * t, Hc, ut=True, F=F)
    return Integrand.imag


def E_3(t):
    Integrand = adia_eng(t, Hc, ut=True, F=F)
    return Integrand.real


def func_psi(t, var):
    return func_psi_module(t, Hc, var)


def Stokes_phase(eps_0):
    delta = adia_param(eps_0, F, D_z, k)

    term1 = math.pi / 4
    term2 = delta * (math.log(delta) - 1)
    term3 = cmath.phase(scipy.special.gamma(1 - 1j * delta))
    return term1 + term2 + term3


for initial_v in [-5, 5]:
    # 基準
    prev_ans_phi_s = Stokes_phase(initial_v)

    if initial_v < 0:  # start from v = -5
        # v_val_find_Stokes is [-2, -4, -6, ...]
        v_val_find_Stokes = eps_0_val[eps_0_val < 0]
        v_val_find_Stokes = v_val_find_Stokes[:: -1]

    else:  # start from v = 5
        # v_val_find_Stokes is [2, 4, 6, ...]
        v_val_find_Stokes = eps_0_val[eps_0_val > 0]

    for eps_0 in v_val_find_Stokes:
        ans_phi_s = 0
        min_error = 10
        for d_phi_s in np.linspace(-0.1, 0.1, 10):
            phi_s = prev_ans_phi_s + d_phi_s  # Stokes phase(弧度法)
            TLZ = -math.pi * (D_z - k*eps_0*F/4)**2 / (abs(eps_0) * abs(F))
            # １回目の遷移がOkaモデルと全体の符号が反転している場合は分子の第２項の符号をマイナスにする
            zero_approx = abs(D_z - k*abs(eps_0)*F/4) / (abs(eps_0) * (-F))

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

            var_init_tmp = eig_vec(t_i, Hc, "upper").tolist()  # initial state
            var_init = [var_init_tmp[0].real, var_init_tmp[0].imag,
                        var_init_tmp[1].real, var_init_tmp[1].imag]
            var_list = solve_ivp(func_psi, [t_i, t_f], var_init, method="LSODA",
                                 t_eval=t_eval, rtol=1e-12, atol=1e-12)
            for i in range(n):
                a = var_list.y[0][i]  # 波動関数 第1成分 実部
                b = var_list.y[1][i]  # 波動関数 第1成分 虚部
                c = var_list.y[2][i]  # 波動関数 第2成分 実部
                d = var_list.y[3][i]  # 波動関数 第2成分 虚部
                psi = np.array([[a + b*1j],
                                [c + d*1j]])  # 波動関数

                # 断熱状態に指定
                q_f = eig_vec(var_list.t[i], Hc, "lower")  # final state

                PA = np.vdot(q_f, psi)  # probability amplitude
                OP = abs(PA)**2  # ocupation probability
                OP_list.append(OP)

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

        print("v = ", eps_0, "Stokes = ", ans_phi_s, "error = ", min_error)
        prev_ans_phi_s = ans_phi_s  # 1つ前の点を基準にする
        Stokes_val.append(ans_phi_s)

    if initial_v < 0:
        Stokes_val = Stokes_val[:: -1]
print("completed")

for eps_0 in eps_0_val:
    Stokes_val_thr_TLZ.append(Stokes_phase(eps_0))

k_tmp = k
k = 0
for eps_0 in eps_0_val:
    Stokes_val_thr_LZ.append(Stokes_phase(eps_0))

k = k_tmp

# %%[markdown]
# グラフの設定

# %%
plt.plot(eps_0_val, Stokes_val, linestyle="None", marker="x", label=rf"$\Delta_y = {"{:.0f}".format(eps_0**2 * k / 4)}$")
plt.plot(eps_0_val, Stokes_val_thr_TLZ, label=rf"$\Delta_y = {"{:.0f}".format(eps_0**2 * k / 4)}$", color="tab:green")
plt.plot(eps_0_val, Stokes_val_thr_LZ, label=r"$\Delta_y = 0$", color="tab:orange")
plt.xlabel(r"energy slope $\varepsilon_0$")
plt.ylabel(r"Stokes phase $\varphi_s$")
plt.title(rf"$\Delta_z = {D_z}, \omega = {-F}$")
plt.legend()
plt.ylim(-0.1, 1.1)
plt.show()

# %%
# Sample
# v_val = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, -5.0, -10.0, -15.0, -20.0, -25.0, -30.0, -35.0, -40.0, -45.0, -50.0, -55.0, -60.0, -65.0, -70.0, -75.0, -80.0, -85.0, -90.0, -95.0, -100.0]
# Stokes_val = [0.20287276377942887, 0.21817888622840845, 0.160015620922286, 0.15083194745289824, 0.23960745765697986, 0.2610360290855513, 0.28858704949371455, 0.31001562092228596, 0.343689090310041, 0.33450541684065327, 0.3314441923508573, 0.4140972535753471, 0.3191992943916736, 0.3957299066365716, 0.36817888622840833, 0.3651176617386124, 0.40491358010595935, 0.4691992943916736, 0.44164827398351036, 0.4508319474528981, 0.20287276377942887, 0.3011486258483944, 0.3477003499863255, 0.44597621205529103, 0.5545969017104635, 0.6218382810208083, 0.6787348327449463, 0.7356313844690843, 0.7821831086070153, 0.8183900051587395, 0.8649417292966706, 0.8908037982621878, 0.9477003499863257, 0.9218382810208084, 0.8752865568828774, 1.0252865568828773, 0.9166658672277048, 0.9942520741242565, 0.9683900051587393, 0.9735624189518428]
