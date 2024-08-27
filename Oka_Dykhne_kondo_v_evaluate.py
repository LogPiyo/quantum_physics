# %%[markdown]
# Oka_Dykhne_kondo.pyで$F$を固定して$\nu$を変化させるプログラムです。
# energy slope $\nu$を横軸、遷移確率$P$を縦軸にしたグラフを出力します。
# $\nu$, $F$の符号反転はそれぞれ、時間反転、エネルギー反転に対応します。
# - 現状Integration Warningが出ます。→調査中

# %%
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad

# parameter
v_values = np.linspace(-2, 2, 100)  # energy slope
m = 0.1  # minimal energy gap
k = 1  # geodesic curvature

# constant
h = 1  # Dirac constant (should not change)
F = -1  # sweep speed (should not change)(時間反転させないため)
TP_list = []  # transition probability


def TLZ_theoretical(v):
    TLZ = -math.pi * (m + k*v*F/4)**2 / (abs(v) * abs(F))
    return np.exp(TLZ)


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

    H['x'] = -v * cmath.cos(q(t))
    H['y'] = 0.25 * k * v**2 * cmath.cos(q(t)) * cmath.sin(2*q(t))
    H['z'] = m * cmath.sin(q(t))
    H['x_dot'] = v * cmath.sin(q(t))
    H['y_dot'] = (0.25 * k * v**2
                  * (-cmath.sin(q(t))*cmath.sin(2*q(t))
                     + 2*cmath.cos(q(t))*cmath.cos(2*q(t))))
    H['z_dot'] = m * cmath.cos(q(t))

    return H[component]


def phi_dot(t):
    """
    define derivative of azimuthal angle with respect to parameter sweep

    Args:
        t (float): time

    Returns:
        float: derivative of azimuthal angle with respect to parameter sweep
    """
    num = -Hc(t, "x")*Hc(t, "y_dot") + Hc(t, "x_dot")*Hc(t, "y")
    den = Hc(t, "x")**2 + Hc(t, "y")**2
    return num / den


def E_plus_unitary_transformed(t):
    """
    define adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: adiabatic energy (unitary transformed)
    """
    X = Hc(t, "x")
    Y = Hc(t, "y")
    Z = Hc(t, "z")
    phi_d = phi_dot(t)
    return cmath.sqrt(X**2 + Y**2 + (Z + 0.5*(-F)*phi_d)**2)


def Re_E(t):
    """
    define real part of adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: real part of adiabatic energy (unitary transformed)
    """
    Integrand = E_plus_unitary_transformed(tp + 1j*t)
    return Integrand.real


for v in v_values:
    tp = math.pi / (2*abs(F))  # transition time
    zero_approx = abs(m - k*abs(v)*F/4) / (abs(v) * (-F))
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

plt.plot(v_values, TP_list, label="numerical")
plt.plot(v_values, TLZ_theoretical(v_values),
         linestyle=":", label="theoretical")
plt.legend()
plt.ylim(-0.1, 1.1)
plt.show()

# %%
