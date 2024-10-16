# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のHamiltonianを用いて，
# Lim, Fuchs and Montambaux(2015)の論文の遷移確率
# $$
# P = \exp \left(-\frac{4}{|w|} \int_0^{\mathrm{Im} t_c} dv
#      \mathrm{Re} (E|_{t=\mathrm{Re} t_c + iv}) \right)
# $$
# を算出したとき，完全トンネルが見られることを確かめるためのプログラムです。
# ただし、ゼロ点のみ遷移点近傍で近似しています。

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
tt = 0  # transition time

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
    define complex Hamiltonian

    Args:
        t (float): time
        component (string): component of vector

    Returns:
        float: specified component
    """
    H = {}

    H['x'] = v * q(t)
    H['y'] = 0.5 * k * v**2 * q(t)**2
    H['z'] = m
    H['x_dot'] = v
    H['y_dot'] = k * v**2 * q(t)
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
        float: real part of adiabatic energy (unitary transformed)
    """
    Integrand = E_plus_unitary_transformed(tt + 1j*t)
    return Integrand.real


for F in F_values:
    zero_approx = abs(m + k*v*F/4) / (abs(v) * (-F))
    # zero of adiabatic energy (approximated)

    # imaginary part of integral of adiabatic energy (unitary transformed)
    ll_Re_E = 0  # lower limit
    ul_Re_E = zero_approx  # upper limit
    log_TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)
    log_TP *= -4 * (-F) / abs(F)

    TP = math.exp(log_TP)  # transition probability
    TP_list.append(TP)

plt.plot(F_values, TP_list, label="numerical")
plt.plot(F_values, TLZ_theoretical(F_values),
         linestyle=":", label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
