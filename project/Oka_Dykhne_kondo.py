# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のTwisted Landau-Zenerモデルが、
# 多重交差するHamiltonianを用いて，
# Lim, Fuchs and Montambaux(2015)の遷移確率
# $$
# P = \exp \left(-\frac{4}{|\delta|} \int_0^{\mathrm{Im} t_c} dv
#      \mathrm{Re} (E|_{t=\mathrm{Re} t_c + iv}) \right)
# $$
# ($\delta:$断熱パラメータ)を算出したとき，完全トンネルが見られることを確かめるためのプログラムです。
# このプログラムでは、掃引速度$F$を横軸、遷移確率$P$を縦軸にしたグラフを出力します。
# ただし、断熱エネルギーのゼロ点の虚部$\mathrm{Im} \, t_c$のみ遷移点近傍で近似した表式を使っています。

# %%
import _pathmagic # noqa
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from my_module.function import TLZ_theoretical, q, adia_eng
from scipy.integrate import quad

# parameter
v = 1  # energy slope
m = 0.1  # minimal energy gap
k = 1  # geodesic curvature

F_values = np.linspace(-2, 2, 100)  # sweep speed

# constant
h = 1  # Dirac constant (should not change)
TP_list = []  # transition probability


def Hc(t, component):
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
    H['y'] = 0.25 * k * v**2 * cmath.cos(q(t, F)) * cmath.sin(2*q(t, F))
    H['z'] = m * cmath.sin(q(t, F))
    H['x_dot'] = v * cmath.sin(q(t, F))
    H['y_dot'] = (0.25 * k * v**2
                  * (-cmath.sin(q(t, F))*cmath.sin(2*q(t, F))
                     + 2*cmath.cos(q(t, F))*cmath.cos(2*q(t, F))))
    H['z_dot'] = m * cmath.cos(q(t, F))

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


for F in F_values:
    tp = math.pi / (2*(-F))  # transition time
    zero_approx = abs(m + k*v*(F)/4) / (abs(v) * (-F))
    # zero of adiabatic energy (approximated)

    # imaginary part of integral of adiabatic energy (unitary transformed)
    ll_Re_E = 0  # lower limit
    ul_Re_E = zero_approx  # upper limit
    log_TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)
    log_TP *= -4 * (-F) / abs(F)

    TP = math.exp(log_TP)
    TP_list.append(TP)

plt.plot(F_values, TP_list, label="numerical")
plt.plot(F_values, TLZ_theoretical(v, F_values, m, k),
         linestyle=":", label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
