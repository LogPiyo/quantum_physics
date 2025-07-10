# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のTwisted Landau-Zenerモデルについて、
# Shr&ouml;dinger方程式の数値微分を行うことで、占有確率の時間発展をプロットする
# プログラムです。初期状態がlower stateのとき、upper stateの占有確率を求めます。

# %%
import _pathmagic # noqa
import math
import numpy as np
import matplotlib.pyplot as plt

from my_module.function import q
from my_module.calculator import calculate_occupation_probability
from numpy.typing import NDArray

# parameter
eps_0 = 1  # energy slope
D_z = 0.1  # minimal energy gap
D_y = 0.25  # twist strength
F = -1  # parameter sweep (should not change)
t_i = -math.pi / abs(F)  # initial time
t_f = math.pi / abs(F)  # final time

tt_1 = -math.pi / (2*abs(F))  # first transition time
tt_2 = math.pi / (2*abs(F))  # second transition time
n = 100  # step
t_eval = np.linspace(t_i, t_f, n)  # time

# constant
h = 1  # Dirac constant (should not change)


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
    H['y'] = 0.5 * (4 * D_y / eps_0**2) * eps_0**2 * q(t, F)**2
    H['z'] = D_z
    H['x_dot'] = eps_0
    H['y_dot'] = (4 * D_y / eps_0**2) * eps_0**2 * q(t, F)
    H['z_dot'] = 0

    return H[component]


OP_list: NDArray = calculate_occupation_probability(H, t_i, t_f, n)

TLZ = -math.pi * (D_z + (4 * D_y / eps_0**2)*eps_0*F/4)**2 / (abs(eps_0) * abs(F))
P_TLZ = math.exp(TLZ) + t_eval*0  # theoretical
arr = np.array(OP_list)
plt.plot(t_eval, arr, label="numerical")
plt.plot(t_eval, P_TLZ, label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
