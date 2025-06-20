# %%[markdown]
# Oka_Dykhne.pyで$F$を固定して$\nu$を変化させるプログラムです。
# energy slope $\nu$を横軸、遷移確率$P$を縦軸にしたグラフを出力します。
# Takayoshi, Wu and Oka(2021)のHamiltonianにおいて、$\nu$を変化させることと
# $F$を変化させることは等価であるため、Oka_Dykhne.pyと見かけ上まったく同じグラフになります。
# ただし解釈はそれぞれ異なります。$\nu$, $F$の符号反転はそれぞれ、時間反転、エネルギー反転に対応します。

# %%
import _pathmagic  # noqa
import math
import numpy as np
import matplotlib.pyplot as plt

from my_module.function import TLZ_theoretical, q, adia_eng, to_LZ
from scipy.integrate import quad

# parameter
v_values = np.linspace(-2, 2, 100)  # energy slope
m = 0.1  # minimal energy gap
k = 1  # geodesic curvature

F = 1  # sweep speed
tt = 0  # transition time

# constant
h = 1  # Dirac constant
TP_list = []  # transition probability


def Hc(t, component):
    """
    define complex Hamiltonian

    Args:
        t (float): time
        component (string): component of vector

    Returns:
        float: specified component
    """
    H = {}

    H['x'] = v * q(t, F)
    H['y'] = 0.5 * k * v**2 * q(t, F)**2
    H['z'] = m
    H['x_dot'] = v
    H['y_dot'] = k * v**2 * q(t, F)
    H['z_dot'] = 0

    return H[component]


def Re_E(t):
    """
    define real part of adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: adiabatic energy
    """
    Integrand = adia_eng(tt + 1j*t, to_LZ(Hc, F))
    return Integrand.real


for v in v_values:
    zero_approx = abs(m + k*v*F/4) / (abs(v) * (-F))
    # zero of adiabatic energy (approximated)

    # imaginary part of integral of adiabatic energy (unitary transformed)
    ll_Re_E = 0  # lower limit
    ul_Re_E = zero_approx  # upper limit
    log_TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)
    log_TP *= -4 * (-F) * abs(F)

    TP = math.exp(log_TP)
    TP_list.append(TP)

plt.plot(v_values, TP_list, label="numerical")
plt.plot(v_values, TLZ_theoretical(v_values, F, m, k),
         linestyle=":", label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
