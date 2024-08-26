# %%[markdown]
# これはTakayoshi, Wu and Oka(2021)のTwisted Landau-Zenerモデルについて、
# Shr&ouml;dinger方程式の数値微分を行うことで、占有確率の時間発展をプロットする
# プログラムです。初期状態がlower stateのとき、upper stateの占有確率を求めます。

# %%
import math
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp

# parameter
v = 1  # energy slope
m = 0.1  # minimal energy gap
k = 1  # geodesic curvature
F = -1  # parameter sweep
t_i = -math.pi / abs(F)  # initial time
t_f = math.pi / abs(F)  # final time

tt_1 = -math.pi / (2*abs(F))  # first transition time
tt_2 = math.pi / (2*abs(F))  # second transition time
n = 200  # step
t_eval = np.linspace(t_i, t_f, n)  # time
OP_list = []  # occupation probability

h = 1  # Dirac constant


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


def eig_vec(t, s):
    """
    eigenvector

    Args:
        t (float): time
        s (state): upper or lower

    Returns:
        array: eigenvector
    """
    E_1 = math.sqrt(H(t, "x")**2 + H(t, "y")**2 + H(t, "z")**2)  # 断熱エネルギー

    # 下の断熱状態を求めるときは断熱エネルギーを符号反転する
    if s == "lower":
        E_1 = -E_1

    eig_vec = np.array([H(t, "x") - H(t, "y")*1j, E_1 - H(t, "z")])
    eig_vec /= np.linalg.norm(eig_vec)  # normalization
    return eig_vec


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
    dadt = (1/h)*(H(t, "x")*var[3] - H(t, "y")*var[2] + H(t, "z")*var[1])
    dbdt = (-1/h)*(H(t, "x")*var[2] + H(t, "y")*var[3] + H(t, "z")*var[0])
    dcdt = (1/h)*(H(t, "x")*var[1] + H(t, "y")*var[0] - H(t, "z")*var[3])
    dddt = (-1/h)*(H(t, "x")*var[0] - H(t, "y")*var[1] - H(t, "z")*var[2])

    return [dadt, dbdt, dcdt, dddt]


# 各Fにおけるpsiの時間発展を計算し，t_fにおけるpsiとFをvar_fに追加する。
var_init_tmp = eig_vec(t_i, "upper").tolist()
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
                    [c+d*1j]])
    q_f = eig_vec(var_list.t[i], "lower")
    dot = np.vdot(q_f, psi)
    OP = abs(dot)**2
    OP_list.append(OP)

TLZ = -math.pi * (m + k*v*F/4)**2 / (v * abs(F))
P_TLZ = math.exp(TLZ) + t_eval*0
arr = np.array(OP_list)
plt.plot(t_eval, arr, label="numerical")
plt.plot(t_eval, P_TLZ, label="theoretical")
plt.ylim(-0.1, 1.1)
