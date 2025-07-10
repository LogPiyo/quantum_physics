# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のTwisted Landau-Zenerモデルについて、
# Shr&ouml;dinger方程式の数値微分を行うことで、占有確率の時間発展をプロットする
# プログラムです。初期状態がlower stateのとき、upper stateの占有確率を求めます。

# %%
# self-made libraries
import _pathmagic # noqa
from my_module.function import q, TLZ_theoretical
from my_module.calculator import calculate_occupation_probability

# type hint libraries
from numpy.typing import NDArray

# math libraries
import math
import numpy as np
import matplotlib.pyplot as plt

# parameter
eps_0 = 1  # energy slope
D_z = 0.1  # minimal energy gap
D_y = 0.25  # twist strength
F = -1  # parameter sweep (should not change)
t_i = -math.pi / abs(F)  # initial time
t_f = math.pi / abs(F)  # final time

n = 100  # step
t_eval = np.linspace(t_i, t_f, n)  # time


def H(t, component):
    hamiltonian = {
        'x': eps_0 * q(t, F),
        'y': 0.5 * (4 * D_y / eps_0**2) * eps_0**2 * q(t, F)**2,
        'z': D_z,
        'x_dot': eps_0,
        'y_dot': (4 * D_y / eps_0**2) * eps_0**2 * q(t, F),
        'z_dot': 0
    }
    return hamiltonian[component]


occupation_prob_arr: NDArray = calculate_occupation_probability(H, t_i, t_f, n)
transition_prob_arr: NDArray = TLZ_theoretical(eps_0, F, D_z, 4 * D_y / eps_0**2) + t_eval * 0

plt.plot(t_eval, occupation_prob_arr, label="numerical")
plt.plot(t_eval, transition_prob_arr, label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
