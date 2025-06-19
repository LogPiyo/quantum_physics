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

import numpy.typing as npt
from typing import Final

# parameter
v: float = 1  # energy slope
m: float = 0.1  # minimal energy gap
k: float = 1  # geodesic curvature
F: float = -1  # parameter sweep (should not change)
t_i: float = -math.pi / abs(F)  # initial time
t_f: float = math.pi / abs(F)  # final time

tt_1: float = -math.pi / (2*abs(F))  # first transition time
tt_2: float = math.pi / (2*abs(F))  # second transition time
n: int = 100  # step
t_eval: npt.NDArray = np.linspace(t_i, t_f, n)  # time

# constant
h: Final[float] = 1  # Dirac constant (should not change)
OP_list: list[float] = []  # occupation probability


def H(t: float, component: str, real: bool = True) -> float:
    """
    define real Hamiltonian

    Args:
        t (float): time
        component (string): component of vector

    Returns:
        float: specified component
    """
    H: dict[str, float] = {}

    H['x'] = v * q(t, F)
    H['y'] = 0.5 * k * v**2 * q(t, F)**2
    H['z'] = m
    H['x_dot'] = v
    H['y_dot'] = k * v**2 * q(t, F)
    H['z_dot'] = 0

    return H[component]


arr: npt.NDArray = calculate_occupation_probability(H, t_i, t_f, n)

TLZ = -math.pi * (m + k*v*F/4)**2 / (abs(v) * abs(F))
P_TLZ = math.exp(TLZ) + t_eval*0  # theoretical
plt.plot(t_eval, arr, label="numerical")
plt.plot(t_eval, P_TLZ, label="theoretical")
plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
