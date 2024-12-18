# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のTwisted Landau-Zenerモデルについて、
# Shr&ouml;dinger方程式の数値微分を行うことで、占有確率の時間発展をプロットする
# プログラムです。初期状態がlower stateのとき、upper stateの占有確率を求めます。

# %%
import _pathmagic # noqa
import math
import numpy as np
import matplotlib.pyplot as plt

from my_module.function import q, eig_vec, func_psi_module
from scipy.integrate import solve_ivp

# parameter
eps_0 = 1  # energy slope
D_z = 0.1  # minimal energy gap
k = 1  # geodesic curvature
F = -1  # parameter sweep (should not change)
t_i = -math.pi / abs(F)  # initial time
t_f = math.pi / abs(F)  # final time

tt_1 = -math.pi / (2*abs(F))  # first transition time
tt_2 = math.pi / (2*abs(F))  # second transition time
n = 100  # step
t_eval = np.linspace(t_i, t_f, n)  # time

# constant
h = 1  # Dirac constant (should not change)
OP_list = []  # occupation probability


def H(t, component, real=True):
    """
    define real Hamiltonian

    Args:
        t (float): time
        component (string): component of vector

    Returns:
        float: specified component
    """
    H = {}

    H['x'] = eps_0 * q(t, F)
    H['y'] = 0.5 * k * eps_0**2 * q(t, F)**2
    H['z'] = D_z
    H['x_dot'] = eps_0
    H['y_dot'] = k * eps_0**2 * q(t, F)
    H['z_dot'] = 0

    return H[component]


def func_psi(t, var):
    return func_psi_module(t, H, var)


# 各Fにおけるpsiの時間発展を計算し，t_fにおけるpsiとFをvar_fに追加する。
var_init_tmp = eig_vec(t_i, H, "upper").tolist()
var_init = [var_init_tmp[0].real, var_init_tmp[0].imag,
            var_init_tmp[1].real, var_init_tmp[1].imag]
var_list = solve_ivp(func_psi, [t_i, t_f], var_init, method="LSODA",
                     t_eval=t_eval, rtol=1e-12, atol=1e-12)

for i in range(n):
    a = var_list.y[0][i]
    b = var_list.y[1][i]
    c = var_list.y[2][i]
    d = var_list.y[3][i]
    psi = np.array([[a+b*1j],
                    [c+d*1j]])  # state vector
    q_f = eig_vec(var_list.t[i], H, "lower")  # final state
    dot = np.vdot(q_f, psi)
    OP = abs(dot)**2  # occupation probability
    OP_list.append(OP)

TLZ = -math.pi * (D_z + k*eps_0*F/4)**2 / (abs(eps_0) * abs(F))
P_TLZ = math.exp(TLZ) + t_eval*0  # theoretical
arr = np.array(OP_list)
plt.plot(t_eval, arr, label="numerical")
plt.plot(t_eval, P_TLZ, label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
