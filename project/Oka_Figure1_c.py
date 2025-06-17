# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のFigure1(c)を再現するプログラムです。
# 任意のHamiltonianで試すことができます。

# %%
import _pathmagic # noqa

import numpy as np
import matplotlib.pyplot as plt

from my_module.function import TLZ_theoretical, q, func_psi_module, eig_vec
from scipy.integrate import solve_ivp

# parameter
v = 1  # energy slope
m = 0.1  # minimal energy gap
k = 1  # geodesic curvature

h = 1  # Dirac constant

n = 50  # step
t_i = -10  # initial time
t_f = 10  # final time

TP_list = []  # transition probability

t_eval = np.linspace(-t_f, t_f, n)  # time
F_values = np.linspace(-2, 2, n)  # sweep speed


def H(t, component):
    """
    define real Hamiltonian

    Args:
        t (float): time
        component (string): 成分

    Returns:
        float: 時刻tにおけるcomponentで指定した成分
    """
    H = {}

    H['x'] = v * q(t, F)
    H['y'] = 0.5 * k * v**2 * q(t, F)**2
    H['z'] = m
    H['x_dot'] = v
    H['y_dot'] = k * v**2 * q(t, F)
    H['z_dot'] = 0

    return H[component]


def func_psi(t, var):
    return func_psi_module(t, H, var)


# 各Fにおけるpsiの時間発展を計算する
for F in F_values:
    var_init_tmp = eig_vec(t_i, H, "upper").tolist()  # initial state
    var_init = [var_init_tmp[0].real, var_init_tmp[0].imag,
                var_init_tmp[1].real, var_init_tmp[1].imag]
    var_list = solve_ivp(func_psi, [t_i, t_f], var_init,
                         method="LSODA", rtol=1e-12, atol=1e-12)
    # 各Fにおける固有関数とpsiを求め，トンネル確率TPを算出する
    a = var_list.y[0][-1]
    b = var_list.y[1][-1]
    c = var_list.y[2][-1]
    d = var_list.y[3][-1]
    psi = np.array([[a+b*1j],
                    [c+d*1j]])  # state vector
    q_f = eig_vec(t_f, H, "lower")  # final state
    dot = np.vdot(q_f, psi)
    TP = abs(dot)**2  # transition probability
    TP_list.append(TP)

plt.plot(F_values, TP_list, label="numerical")
plt.plot(F_values, TLZ_theoretical(v, F_values, m, k),
         linestyle=":", label="theoretical")
plt.legend()
plt.show()
