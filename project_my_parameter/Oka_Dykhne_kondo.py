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

from my_module.function import TLZ_theoretical, q, adia_eng, to_LZ
from scipy.integrate import quad

# parameter
eps_0 = 1  # energy slope
D_z = 0.1  # minimal energy gap
D_y = 0.25  # twist strength

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

    H['x'] = -eps_0 * cmath.cos(q(t, F))
    H['y'] = 0.25 * (4 * D_y / eps_0**2) * eps_0**2 * cmath.cos(q(t, F)) * cmath.sin(2*q(t, F))
    H['z'] = D_z * cmath.sin(q(t, F))
    H['x_dot'] = eps_0 * cmath.sin(q(t, F))
    H['y_dot'] = (0.25 * (4 * D_y / eps_0**2) * eps_0**2
                  * (-cmath.sin(q(t, F))*cmath.sin(2*q(t, F))
                     + 2*cmath.cos(q(t, F))*cmath.cos(2*q(t, F))))
    H['z_dot'] = D_z * cmath.cos(q(t, F))

    return H[component]


def Re_E(t):
    """
    define real part of adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: real part of adiabatic energy (unitary transformed)
    """
    Integrand = adia_eng(tp + 1j*t, to_LZ(Hc, F))
    return Integrand.real


for F in F_values:
    tp = math.pi / (2*(-F))  # transition time
    zero_approx = abs(D_z + (4 * D_y / eps_0**2)*eps_0*(F)/4) / (abs(eps_0) * (-F))
    # zero of adiabatic energy (approximated)

    # imaginary part of integral of adiabatic energy (unitary transformed)
    ll_Re_E = 0  # lower limit
    ul_Re_E = zero_approx  # upper limit
    log_TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)
    log_TP *= -4 * (-F) / abs(F)

    TP = math.exp(log_TP)
    TP_list.append(TP)

plt.plot(F_values, TP_list, label="numerical")
plt.plot(F_values, TLZ_theoretical(eps_0, F_values, D_z, (4 * D_y / eps_0**2)),
         linestyle=":", label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
