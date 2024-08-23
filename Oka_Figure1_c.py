# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のFigure1(c)を再現するプログラムです。
# 任意のHamiltonianで試すことができます。

# %%
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

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

    H['x'] = v * q(t)
    H['y'] = 0.5 * k * v**2 * q(t)**2
    H['z'] = m
    H['x_dot'] = v
    H['y_dot'] = k * v**2 * q(t)
    H['z_dot'] = 0

    return H[component]


def E_adia(t):
    """
    adiabatic energy (eigenvalue)

    Args:
        t (float): time

    Returns:
        float: adiabatic enegy
    """
    E_adia = math.sqrt(H(t, "x")**2 + H(t, "y")**2 + H(t, "z")**2)
    return E_adia


def func_psi(t, var):
    """
    state vector

    (t_f)における系の状態ベクトル(psi(t_f))を求める関数です。
    psiの第1成分をa+ib, 第2成分をc*idとします。
    var[0]=a,var[1]=b, var[2]=c, var[3]=dとします。

    Args:
        t (float): time
        var (list): 状態ベクトルの各成分を要素とするlist

    Returns:
        list: 微分方程式
    """
    dadt = (1/h) * (H(t, "x")*var[3] - H(t, "y")*var[2] + H(t, "z")*var[1])
    dbdt = (-1/h) * (H(t, "x")*var[2] + H(t, "y")*var[3] + H(t, "z")*var[0])
    dcdt = (1/h) * (H(t, "x")*var[1] + H(t, "y")*var[0] - H(t, "z")*var[3])
    dddt = (-1/h) * (H(t, "x")*var[0] - H(t, "y")*var[1] - H(t, "z")*var[2])

    return [dadt, dbdt, dcdt, dddt]


def eig_vec(t, s):
    """
    eigenvector

    Args:
        t (float): time
        s (string): select upper or lower state

    Returns:
        array: eigenvector
    """
    energy = E_adia(t)
    # 下の断熱状態を求めるときは断熱エネルギーを符号反転する
    if s == "lower":
        energy = -energy

    eig_vec = np.array([H(t, "x") - H(t, "y")*1j, energy - H(t, "z")])
    eig_vec /= np.linalg.norm(eig_vec)  # normalization
    return eig_vec


# 各Fにおけるpsiの時間発展を計算する
for F in F_values:
    var_init_tmp = eig_vec(t_i, "upper").tolist()  # initial state
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
    q_f = eig_vec(t_f, "lower")  # final state
    dot = np.vdot(q_f, psi)
    TP = abs(dot)**2  # transition probability
    TP_list.append(TP)

plt.plot(F_values, TP_list)
plt.show()

# %%
