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
from scipy.optimize import newton

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
F_eval = np.linspace(-2, 2, n)  # time


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
    H['x_dot_dot'] = 0
    H['y_dot_dot'] = k * v**2
    H['z_dot_dot'] = 0

    return H[component]


def E_plus(t):
    """
    adiabatic energy (eigenvalue)

    Args:
        t (float): time

    Returns:
        float: adiabatic enegy
    """
    E_plus = cmath.sqrt(Hc(t, "x")**2 + Hc(t, "y")**2 + Hc(t, "z")**2)
    return E_plus


def phi_dot(t):
    """
    define first derivative of phi

    Args:
        t (float): time

    Returns:
        float: first derivative of phi
    """
    num = -Hc(t, "x")*Hc(t, "y_dot") + Hc(t, "x_dot")*Hc(t, "y")
    den = Hc(t, "x")**2 + Hc(t, "y")**2
    return num / den


def phi_dot_dot(t):
    """
    define second derivative of phi

    Args:
        t (float): time

    Returns:
        float: second derivative of phi
    """
    X = Hc(tt + 1j*t, "x")
    Y = Hc(tt + 1j*t, "y")
    X_d = Hc(tt + 1j*t, "x_dot")
    Y_d = Hc(tt + 1j*t, "y_dot")
    X_dd = Hc(tt + 1j*t, "x_dot_dot")
    Y_dd = Hc(tt + 1j*t, "y_dot_dot")

    num_term1 = (X*Y_dd - Y*X_dd) * (X**2 + Y**2)
    num_term2 = 2 * (X*Y_d - Y*X_d) * (X*X_d + Y*Y_d)
    num = num_term1 - num_term2
    den = (X**2 + Y**2)**2

    return -num / den


def Re_E(t):
    """
    define real part of adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: adiabatic energy
    """
    X = Hc(tt + 1j*t, "x")
    Y = Hc(tt + 1j*t, "y")
    Z = Hc(tt + 1j*t, "z")
    phi_d = phi_dot(tt + 1j*t)
    a = X**2 + Y**2 + (Z + 0.5*(-F)*phi_d)**2
    Integrand = (cmath.sqrt(a))
    return -4 * (-F) / abs(F) * Integrand.real


def E_for_zero(t):
    X = Hc(tt + 1j*t, "x")
    Y = Hc(tt + 1j*t, "y")
    Z = Hc(tt + 1j*t, "z")
    phi_d = phi_dot(tt + 1j*t)

    return X**2 + Y**2 + (Z + 0.5*(-F)*phi_d)**2


def derv_E_for_zero(t):
    """
    derivative of adiabatic energy (unitary transformed) for finding zero

    Note: 正確にはこの関数の戻り値の平方根が断熱エネルギーである。

    Args:
        t (float): time

    Returns:
        float: derivative of adiabatic energy
    """
    X = Hc(tt + 1j*t, "x")
    Y = Hc(tt + 1j*t, "y")
    Z = Hc(tt + 1j*t, "z")
    X_dot = Hc(tt + 1j*t, "x_dot")
    Y_dot = Hc(tt + 1j*t, "y_dot")
    Z_dot = Hc(tt + 1j*t, "z_dot")
    phi_d = phi_dot(tt + 1j*t)
    phi_dd = phi_dot_dot(tt + 1j*t)

    term1 = 2 * X * X_dot
    term2 = 2 * Y * Y_dot
    term3 = 2 * (Z - 0.5*(-F)*phi_d) * (Z_dot - 0.5*(-F)*phi_dd)

    return term1 + term2 + term3


for F in F_eval:
    tt = 0.01  # transition time
    zero_approx = (m + k*v*abs(F)/4) / (v * (-F))
    closest_zero = newton(E_for_zero, tt+1j*zero_approx, fprime=derv_E_for_zero, maxiter=1000)
    print(closest_zero)
    zero_approx = (m + k*v*abs(F)/4) / (v * (-F))
    # zero of adiabatic energy (approximated)

    # imaginary part of integral of adiabatic energy (unitary transformed)
    ll_Re_E = 0  # lower limit
    ul_Re_E = zero_approx  # upper limit
    log_TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)

    TP = math.exp(log_TP)
    TP_list.append(TP)

plt.plot(F_eval, TP_list)
plt.ylim(-0.1, 1.1)
plt.show()

# %%
