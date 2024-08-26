# %%[markdown]のHamiltonianを用いて，
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

t_i = -math.pi  # initial time
t_f = math.pi  # final time

# constant
h = 1  # Dirac constant
TP_list = []  # transition probability
n = 100  # step
F_values = np.linspace(-2, 2, n)  # time


def TLZ_theoretical(F):
    TLZ = -math.pi * (m + k*v*F/4)**2 / (v * abs(F))
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
        component (string): 取り出したい成分

    Returns:
        float: 時刻tにおけるcomponentで指定した成分
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
    define derivative of phi

    Args:
        t (float): time

    Returns:
        float: derivative of phi
    """
    num = -Hc(t, "x")*Hc(t, "y_dot") + Hc(t, "x_dot")*Hc(t, "y")
    den = Hc(t, "x")**2 + Hc(t, "y")**2
    return num / den


def E_plus_unitary_transformed(t):
    """
    adiabatic energy (unitary transformed)

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
    return -4 * (-F) / abs(F) * Integrand.real


for F in F_values:
    tt = 0  # transition time
    zero_approx = (m + k*v*abs(F)/4) / (v * (-F))
    # zero of adiabatic energy (approximated)

    # imaginary part of integral of adiabatic energy (unitary transformed)
    ll_Re_E = 0  # lower limit
    ul_Re_E = zero_approx  # upper limit
    log_TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)

    TP = math.exp(log_TP)
    TP_list.append(TP)

plt.plot(F_values, TP_list, label="numerical")
plt.plot(F_values, TLZ_theoretical(F_values),
         linestyle=":", label="theoretical")
plt.legend()
plt.ylim(-0.1, 1.1)
plt.show()
