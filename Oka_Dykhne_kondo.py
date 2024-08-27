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
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad

# parameter
v = 1  # energy slope
m = 0.1  # minimal energy gap
k = 1  # geodesic curvature

F_values = np.linspace(-2, 2, 100)  # sweep speed

# constant
h = 1  # Dirac constant (should not change)
TP_list = []  # transition probability


def TLZ_theoretical(F):
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
plt.plot(F_values, TLZ_theoretical(F_values),
         linestyle=":", label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
