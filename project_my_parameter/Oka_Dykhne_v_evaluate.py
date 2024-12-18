# %%[markdown]
# Oka_Dykhne.pyで$F$を固定して$\nu$を変化させるプログラムです。
# energy slope $\nu$を横軸、遷移確率$P$を縦軸にしたグラフを出力します。
# Takayoshi, Wu and Oka(2021)のHamiltonianにおいて、$\nu$を変化させることと
# $F$を変化させることは等価であるため、Oka_Dykhne.pyと見かけ上まったく同じグラフになります。
# ただし解釈はそれぞれ異なります。$\nu$, $F$の符号反転はそれぞれ、時間反転、エネルギー反転に対応します。

# %%
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad

# parameter
eps_0_values = np.linspace(-2, 2, 100)  # energy slope
D_z = 0.1  # minimal energy gap
k = 1  # geodesic curvature

F = 1  # sweep speed
tt = 0  # transition time

# constant
h = 1  # Dirac constant
TP_list = []  # transition probability


def TLZ_theoretical(v):
    TLZ = -math.pi * (D_z + k*v*F/4)**2 / (abs(v) * abs(F))
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
    define complex Hamiltonian

    Args:
        t (float): time
        component (string): component of vector

    Returns:
        float: specified component
    """
    H = {}

    H['x'] = eps_0 * q(t)
    H['y'] = 0.5 * k * eps_0**2 * q(t)**2
    H['z'] = D_z
    H['x_dot'] = eps_0
    H['y_dot'] = k * eps_0**2 * q(t)
    H['z_dot'] = 0

    return H[component]


def phi_dot(t):
    """
    define derivative of azimuthal angle

    Args:
        t (float): time

    Returns:
        float: derivative of azimuthal angle
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
        float: adiabatic energy
    """
    Integrand = E_plus_unitary_transformed(tt + 1j*t)
    return Integrand.real


for eps_0 in eps_0_values:
    zero_approx = abs(D_z + k*eps_0*F/4) / (abs(eps_0) * (-F))
    # zero of adiabatic energy (approximated)

    # imaginary part of integral of adiabatic energy (unitary transformed)
    ll_Re_E = 0  # lower limit
    ul_Re_E = zero_approx  # upper limit
    log_TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)
    log_TP *= -4 * (-F) * abs(F)

    TP = math.exp(log_TP)
    TP_list.append(TP)

plt.plot(eps_0_values, TP_list, label="numerical")
plt.plot(eps_0_values, TLZ_theoretical(eps_0_values),
         linestyle=":", label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
