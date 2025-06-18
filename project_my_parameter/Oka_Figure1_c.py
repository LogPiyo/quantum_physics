# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のFigure1(c)を再現するプログラムです。
# 任意のHamiltonianで試すことができます。

# %%
import _pathmagic # noqa

from typing import Final
import numpy.typing as npt

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from my_module.function import TLZ_theoretical, q, func_psi_module, eig_vec

# parameter
eps_0: float = 1  # energy slope
D_z: float = 0.1  # minimal energy gap
D_y: float = 0.25  # twist strength

h: Final[float] = 1  # Dirac constant

n: int = 50  # step
t_i: float = -10  # initial time
t_f: float = 10  # final time

TP_list: list[float] = []  # transition probability

t_eval: npt.NDArray = np.linspace(-t_f, t_f, n)  # time
F_values: npt.NDArray = np.linspace(-2, 2, n)  # sweep speed


def H(t: float, component: str) -> float:
    """
    define real Hamiltonian

    Args:
        t (float): time
        component (string): 成分

    Returns:
        float: 時刻tにおけるcomponentで指定した成分
    """
    H: dict[str, float] = {}

    H['x'] = eps_0 * q(t, F)
    H['y'] = 0.5 * (4 * D_y / eps_0**2) * eps_0**2 * q(t, F)**2
    H['z'] = D_z
    H['x_dot'] = eps_0
    H['y_dot'] = (4 * D_y / eps_0**2) * eps_0**2 * q(t, F)
    H['z_dot'] = 0

    return H[component]


def func_psi(t: float, var: list[float]) -> list[float]:
    return func_psi_module(t, H, var)


# 各Fにおけるpsiの時間発展を計算する
for F in F_values:
    var_init_tmp: list[float] = eig_vec(t_i, H, "upper").tolist()  # initial state
    var_init: list[float] = [var_init_tmp[0].real, var_init_tmp[0].imag,
                             var_init_tmp[1].real, var_init_tmp[1].imag]
    var_list = solve_ivp(func_psi, [t_i, t_f], var_init,
                         method="LSODA", rtol=1e-12, atol=1e-12)  # 型はscipy.integrate.OdeResult
    # 各Fにおける固有関数とpsiを求め，トンネル確率TPを算出する
    a: float = var_list.y[0][-1]
    b: float = var_list.y[1][-1]
    c: float = var_list.y[2][-1]
    d: float = var_list.y[3][-1]
    psi: npt.NDArray = np.array([[a+b*1j],
                                 [c+d*1j]])  # state vector
    q_f: npt.NDArray = eig_vec(t_f, H, "lower")  # final state
    dot: float = np.vdot(q_f, psi)
    TP: float = abs(dot)**2  # transition probability
    TP_list.append(TP)

# %% グラフの描画
plt.plot(F_values, TP_list, label="numerical")
plt.plot(F_values, TLZ_theoretical(eps_0, F_values, D_z, (4 * D_y / eps_0**2)),
         linestyle=":", label="theoretical")
plt.legend()
plt.show()
