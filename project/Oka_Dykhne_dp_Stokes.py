# %%[markdown]
# これは解析的に求めた占有確率において，Stokes位相をフィッティングパラメータとしたときの
# Stokes位相のエネルギー依存性を可視化するプログラムです。

# %%
# multiple-passage TLZ modelにおける断熱状態の占有確率を数値計算結果と理論値で比較する
import math
import cmath
import scipy
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp, quad


# parameter
v_val = np.linspace(-100, 100, 50)  # energy slope
v_val = v_val[(v_val >= 5) | (v_val <= -5)]
m = 4  # minimal energy gap
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
Stokes_val_LZ = []
t_eval = np.linspace(t_i, t_f, n)  # time


def q(t):
    """
    define parameter sweep q

    q = adiabatic parameter * time

    Args:
        t (float): time

    Returns:
        float: q
    """
    return -F * t


def H(t, component):
    """Hamiltonianの設定(実数のみ)

    Args:
        t (float): time
        component (string): 成分

    Returns:
        float: 時刻tにおけるcomponentで指定した成分を返す。
    """
    H = {}

    H['x'] = -v * math.cos(q(t))
    H['y'] = -0.125 * k * v**2 * math.sin(2 * q(t))**2
    H['z'] = m * math.sin(q(t))

    return H[component]


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

    H['x'] = -v * cmath.cos(q(t))
    H['y'] = -0.125 * k * v**2 * cmath.sin(2 * q(t))**2
    H['z'] = m * cmath.sin(q(t))
    H['x_dot'] = v * cmath.sin(q(t))
    H['y_dot'] = -0.125 * k * v**2 * 4 * cmath.sin(2 * q(t)) * cmath.cos(2 * q(t))
    H['z_dot'] = m * cmath.cos(q(t))

    return H[component]


def E_plus(t):
    """
    adiabatic energy

    Args:
        t (float): time

    Returns:
        float: adiabatic enegy
    """
    E_plus = cmath.sqrt(Hc(t, "x")**2 + Hc(t, "y")**2 + Hc(t, "z")**2)
    return E_plus


def phi_dot(t):
    """
    approximate form of phi dot

    Args:
        t (float): time

    Returns:
        float: phi dot (approximated)
    """
    num = -Hc(t, "x") * Hc(t, "y_dot") + Hc(t, "x_dot") * Hc(t, "y")
    den = Hc(t, "x")**2 + Hc(t, "y")**2
    return num / den


def E_plus_unitary_transformed(t):
    """
    adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: adiabatic enegy (unitary transformed)
    """
    X = Hc(t, "x")
    Y = Hc(t, "y")
    Z = Hc(t, "z")
    phi_d = phi_dot(t)
    return cmath.sqrt(X**2 + Y**2 + (Z + 0.5 * (-F) * phi_d)**2)


def Re_E(t):
    """
    define real part of adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: adiabatic energy
    """
    Integrand = E_plus_unitary_transformed(tp_1 + 1j * t)
    return Integrand.real


def Im_E_1(t):
    Integrand = E_plus_unitary_transformed(tp_1 + 1j * t)
    return Integrand.imag


def Im_E_2(t):
    Integrand = E_plus_unitary_transformed(tp_2 + 1j * t)
    return Integrand.imag


def E_3(t):
    Integrand = E_plus_unitary_transformed(t)
    return Integrand.real


def eig_vec(t, s):
    """
    eigenvector

    Args:
        t (float): time
        s (state): upper or lower

    Returns:
        array: eigenvector
    """
    energy = math.sqrt(H(t, "x")**2 + H(t, "y")**2 + H(t, "z")**2)  # 断熱エネルギー

    # 下の断熱状態を求めるときは断熱エネルギーを符号反転する
    if s == "lower":
        energy = -energy

    eig_vec = np.array([H(t, "x") - H(t, "y") * 1j, energy - H(t, "z")])
    eig_vec /= np.linalg.norm(eig_vec)  # normalization
    return eig_vec


def func_psi(t, var):
    """
    state vector

    (t_f)における系の状態ベクトル(psi(t_f))を求める関数です。
    psiの第1成分をa+ib, 第2成分をc*idとします。
    var[0]=a,var[1]=b, var[2]=c, var[3]=dとします。

    Args:
        t (float): time
        var (list): 状態ベクトルの各成分を要素とするlist

    Returns:
        list: 微分方程式
    """
    dadt = (1/h)*(H(t, "x")*var[3] - H(t, "y")*var[2] + H(t, "z")*var[1])
    dbdt = (-1/h)*(H(t, "x")*var[2] + H(t, "y")*var[3] + H(t, "z")*var[0])
    dcdt = (1/h)*(H(t, "x")*var[1] + H(t, "y")*var[0] - H(t, "z")*var[3])
    dddt = (-1/h)*(H(t, "x")*var[0] - H(t, "y")*var[1] - H(t, "z")*var[2])

    return [dadt, dbdt, dcdt, dddt]


def adia_param(v):
    return m**2 / (2 * abs(v) * abs(F))


def Stokes_phase(v):
    delta = adia_param(v)
    term1 = math.pi / 4
    term2 = delta * (math.log(delta) - 1)
    term3 = cmath.phase(scipy.special.gamma(1 - 1j * delta))
    return term1 + term2 + term3


for initial_v in [-5, 5]:
    # 基準
    prev_ans_phi_s = Stokes_phase(initial_v)

    if initial_v < 0:  # start from v = -5
        # v_val_find_Stokes is [-2, -4, -6, ...]
        v_val_find_Stokes = v_val[v_val < 0]
        v_val_find_Stokes = v_val_find_Stokes[:: -1]

    else:  # start from v = 5
        # v_val_find_Stokes is [2, 4, 6, ...]
        v_val_find_Stokes = v_val[v_val > 0]

    for v in v_val_find_Stokes:
        ans_phi_s = 0
        min_error = 10
        for d_phi_s in np.linspace(-0.1, 0.1, 10):
            phi_s = prev_ans_phi_s + d_phi_s  # Stokes phase(弧度法)
            TLZ = -math.pi * (m - k*v*F/4)**2 / (abs(v) * abs(F))
            # １回目の遷移がOkaモデルと全体の符号が反転している場合は分子の第２項の符号をマイナスにする
            zero_approx = abs(m - k*abs(v)*F/4) / (abs(v) * (-F))

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

            var_init_tmp = eig_vec(t_i, "upper").tolist()  # initial state
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
                q_f = eig_vec(var_list.t[i], "lower")  # final state

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

        print("v = ", v, "Stokes = ", ans_phi_s, "error = ", min_error)
        prev_ans_phi_s = ans_phi_s  # 1つ前の点を基準にする
        Stokes_val.append(ans_phi_s)

    if initial_v < 0:
        Stokes_val = Stokes_val[:: -1]
print("completed")

for v in v_val:
    Stokes_val_LZ.append(Stokes_phase(v))

# %%[markdown]
# グラフの設定

# %%
plt.plot(v_val, Stokes_val, linestyle="None", marker="x", label=r"$D_y = 45$")
plt.plot(v_val, Stokes_val_LZ, label=r"$D_y = 0$")
plt.xlabel(r"energy slope $\epsilon_0$")
plt.ylabel(r"Stokes phase $\psi_s$")
plt.legend()
plt.ylim(-0.1, 1.1)
plt.show()

# %%
# Sample
# v_val = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, -5.0, -10.0, -15.0, -20.0, -25.0, -30.0, -35.0, -40.0, -45.0, -50.0, -55.0, -60.0, -65.0, -70.0, -75.0, -80.0, -85.0, -90.0, -95.0, -100.0]
# Stokes_val = [0.20287276377942887, 0.21817888622840845, 0.160015620922286, 0.15083194745289824, 0.23960745765697986, 0.2610360290855513, 0.28858704949371455, 0.31001562092228596, 0.343689090310041, 0.33450541684065327, 0.3314441923508573, 0.4140972535753471, 0.3191992943916736, 0.3957299066365716, 0.36817888622840833, 0.3651176617386124, 0.40491358010595935, 0.4691992943916736, 0.44164827398351036, 0.4508319474528981, 0.20287276377942887, 0.3011486258483944, 0.3477003499863255, 0.44597621205529103, 0.5545969017104635, 0.6218382810208083, 0.6787348327449463, 0.7356313844690843, 0.7821831086070153, 0.8183900051587395, 0.8649417292966706, 0.8908037982621878, 0.9477003499863257, 0.9218382810208084, 0.8752865568828774, 1.0252865568828773, 0.9166658672277048, 0.9942520741242565, 0.9683900051587393, 0.9735624189518428]
