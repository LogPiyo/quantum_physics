# %%[markdown]
# 遷移点付近のHamiltonianがすべて同じになる系の場合で下記を実現します。
#
# Oka_Dykhne_kondo.pyで$F$を固定して$\nu$を変化させるプログラムです。
# energy slope $\nu$を横軸、遷移確率$P$を縦軸にしたグラフを出力します。
# $\nu$, $F$の符号反転はそれぞれ、時間反転、エネルギー反転に対応します。
# - 現状Integration Warningが出ます。→調査中

# %%
import _pathmagic # noqa
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from my_module.function import q, adia_eng
from scipy.integrate import quad

# parameter
v_values = np.linspace(-200, 200, 50)  # energy slope
m = 4  # minimal energy gap
k = 0.1  # geodesic curvature

# constant
h = 1  # Dirac constant (should not change)
F = -1  # sweep speed (should not change)(時間反転させないため)
TP_list = []  # transition probability


def TLZ_theoretical(v):
    TLZ = -math.pi * (m - k * v * F / 4)**2 / (abs(v) * abs(F))
    return np.exp(TLZ)


def Hc(t, component, real=False):
    """
    define complex Hamiltonian and its derivative
    with respect to parameter sweep

    Args:
        t (float): time
        component (string): component of vector

    Returns:
        float: specified component
    """
    H = {}

    H['x'] = -v * cmath.cos(q(t, F))
    # H['y'] = (0.25 * k * v**2 * cmath.cos(q(t, F)) * cmath.sin(2 * q(t, F)))
    H['y'] = -0.125 * k * v**2 * cmath.sin(2 * q(t, F))*cmath.sin(2 * q(t, F))
    H['z'] = m * cmath.sin(q(t, F))
    H['x_dot'] = v * cmath.sin(q(t, F))
    # H['y_dot'] = (0.25 * k * v**2
    #               * (-cmath.sin(q(t, F)) * cmath.sin(2 * q(t, F))
    #                  + 2*cmath.cos(q(t, F)) * cmath.cos(2 * q(t, F))))
    H['y_dot'] = -0.125 * k * v**2 * 4 * cmath.sin(2 * q(t, F)) * cmath.cos(2 * q(t, F))
    H['z_dot'] = m * cmath.cos(q(t, F))

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
        float: real part of adiabatic energy (unitary transformed)
    """
    Integrand = adia_eng(tp + 1j*t, Hc, ut=True, F=F)
    return Integrand.real


for v in v_values:
    # if abs(v) < 0.07:  # avoid Integration Warning
    #     continue

    tp = -math.pi / (2 * abs(F))  # transition time
    zero_approx = abs(m - k * abs(v) * F / 4) / (abs(v) * (-F))
    # -pi/2のときは分子の第二項の符号が変わる
    # 被積分関数の符号と合わせる
    # zero of adiabatic energy (approximated)

    # imaginary part of integral of adiabatic energy (unitary transformed)
    ll_Re_E = 0  # lower limit
    ul_Re_E = zero_approx  # upper limit
    log_TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)
    log_TP *= -4 * (-F) / abs(F)

    TP = math.exp(log_TP)
    TP_list.append(TP)

# %%
# v_values = v_values[abs(v_values) >= 0.07]  # avoid Integration Warning
plt.plot(v_values, TP_list, label="numerical", color="tab:blue", linestyle="None", marker="x")
plt.plot(v_values, TLZ_theoretical(v_values), label="theoretical", color="tab:red")
plt.xlabel(r"energy slope $\epsilon_0$")
plt.ylabel(r"transition probability $P$")
plt.title(rf"$\Delta_Z = {m}, \omega = {-F}$")
plt.legend()
plt.ylim(-0.1, 1.1)
plt.show()

# %%
