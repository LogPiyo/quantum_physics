# %%[markdown]
# Oka(2021)をもとに
# ユニタリ変換後のHamiltonianで占有確率を数値計算<br>
# $$\Delta_x \sin{\omega t} \sigma_x
#  + Delta_y \cos{\omega t} \sin{2 \omega t}\sigma_y
#  + \epsilon_0 \cos{\omega t} \sigma_z\\$$
# に対してユニタリ変換する

# %%
# multiple-passage TLZ modelにおける断熱状態の占有確率を数値計算結果と理論値で比較する
import math
import cmath
import scipy
from scipy.integrate import solve_ivp, quad
import numpy as np
import matplotlib.pyplot as plt

# parameter
v = 1  # energy slope
m = 0.1  # minimal energy gap
k = 1  # geodesic curvature
F = 1
t_i = -math.pi  # initial time
t_f = math.pi  # final time

tp_1 = -math.pi / (2*(-F))
tp_2 = math.pi / (2*(-F))

# constant
h = 1  # Dirac constant
TP_list = []  # transition probability
n = 100  # step
F_eval = np.linspace(-2, 2, n)  # time

OP_list = []  # ocupation probability
t_eval = np.linspace(t_i, t_f, n)  # time
delta = m**2 / (2 * v * abs(F))  # adiabatic parameter
phi_s = (math.pi/4
         + delta*(math.log(delta)-1)
         + cmath.phase(scipy.special.gamma(1-1j*delta)))  # Stokes phase(弧度法)

TLZ = -math.pi * (m + k*v*F/4)**2 / (v * abs(F))
zero_approx = (m + k*v*abs(F)/4) / (v * (-F))


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
    """Hamiltonianの設定(実数のみ)

    Args:
        t (float): time
        component (string): 成分

    Returns:
        float: 時刻tにおけるcomponentで指定した成分を返す。
    """
    H = {}

    H['x'] = -v * math.cos(q(t))
    H['y'] = 0.25 * k * v**2 * math.cos(q(t)) * math.sin(2*q(t))
    H['z'] = m * math.sin(q(t))
    H['x_dot'] = v * math.sin(q(t))
    H['y_dot'] = (0.25 * k * v**2
                  * (-math.sin(q(t))*math.sin(2*q(t))
                     + 2*math.cos(q(t))*math.cos(2*q(t))))
    H['z_dot'] = m * math.cos(q(t))

    return H[component]


def Hc(t, component):
    """
    Hamiltonianの設定(複素数対応)

    Args:
        t (float): time
        component (string): 成分

    Returns:
        float: 時刻tにおけるcomponentで指定した成分を返す。
    """
    H = {}

    H['x'] = -v * cmath.cos(q(t))
    H['y'] = 0.25 * k * v**2 * cmath.cos(q(t)) * cmath.sin(2*q(t))
    H['z'] = m * cmath.sin(q(t))
    H['x_dot'] = v * cmath.sin(q(t))
    H['y_dot'] = (0.25 * k * v**2
                  * (-cmath.sin(q(t))*cmath.sin(2*q(t))
                     + 2*cmath.cos(q(t))*cmath.cos(2*q(t))))
    H['z_dot'] = m * cmath.cos(q(t))

    return H[component]


def E_plus(t):
    """
    adiabatic energy

    Args:
        t (float): time

    Returns:
        float: adiabatic enegy
    """
    E_plus = cmath.sqrt(Hc(t, "x")**2 + Hc(t, "y")**2 + Hc(t, "z")**2)
    return E_plus


def phi_dot(t):
    """
    approximate form of phi dot

    Args:
        t (float): time

    Returns:
        float: phi dot (approximated)
    """
    num = -Hc(t, "x")*Hc(t, "y_dot") + Hc(t, "x_dot")*Hc(t, "y")
    den = Hc(t, "x")**2 + Hc(t, "y")**2
    return num / den
#     return -k * v / 2


def E_plus_unitary_transformed(t):
    """
    adiabatic energy (unitary transformed)

    Args:
        t (float): time

    Returns:
        float: adiabatic enegy (unitary transformed)
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
    Integrand = E_plus_unitary_transformed(tp_2 + 1j*t)
    return Integrand.real


# integral
ll_Re_E = 0  # lower limit
ul_Re_E = zero_approx  # upper limit
TP, _ = quad(Re_E, ll_Re_E, ul_Re_E)
TP *= -4 * (-F) / abs(F)


def Im_E_1(t):
    Integrand = E_plus_unitary_transformed(tp_1 + 1j*t)
    return Integrand.imag


# integral
ll_Im_E_1 = 0  # lower limit
ul_Im_E_1 = tp_1  # upper limit
phase_term1, err = quad(Im_E_1, ll_Im_E_1, ul_Im_E_1)
phase_term1 *= (-F) / abs(F)


def Im_E_2(t):
    Integrand = E_plus_unitary_transformed(tp_2 + 1j*t)
    return Integrand.imag


# integral
ll_Im_E_2 = 0  # lower limit
ul_Im_E_2 = tp_2  # upper limit
phase_term2, err = quad(Im_E_2, ll_Im_E_2, ul_Im_E_2)
phase_term2 *= (-F) / abs(F)


def E_3(t):
    Integrand = E_plus_unitary_transformed(t)
    return Integrand.real


# integral
ll_E_3 = tp_1  # lower limit
ul_E_3 = tp_2  # upper limit
phase_term3, err = quad(E_3, ll_E_3, ul_E_3)
phase_term3 *= (-F) / abs(F)


def eig_vec(t, s):
    """
    eigenvector

    Args:
        t (float): time
        s (state): plus or minus

    Returns:
        array: eigenvector
    """
    E_1 = math.sqrt(H(t, "x")**2 + H(t, "y")**2 + H(t, "z")**2)  # 断熱エネルギー

    # 下の断熱状態を求めるときは断熱エネルギーを符号反転する
    if s == "minus":
        E_1 = -E_1

    eig_vec = np.array([H(t, "x") - H(t, "y")*1j, E_1 - H(t, "z")])
    eig_vec /= np.linalg.norm(eig_vec)  # normalization
    return eig_vec


def func_psi(t, var):
    """
    state vector

    (t_f)における系の状態ベクトル(psi(t_f))を求める関数です。
    psiの第1成分をa+ib，第2成分をc*idとします。
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


# 各時間における波動関数を算出
var_init_tmp = eig_vec(t_i, "minus").tolist()  # initial state
var_init = [var_init_tmp[0].real, var_init_tmp[0].imag,
            var_init_tmp[1].real, var_init_tmp[1].imag]
var_list = solve_ivp(func_psi, [t_i, t_f], var_init, method="LSODA",
                     t_eval=t_eval, rtol=1e-12, atol=1e-12)
for i in range(n):
    a = var_list.y[0][i]  # 波動関数 第1成分 実部
    b = var_list.y[1][i]  # 波動関数 第1成分 虚部
    c = var_list.y[2][i]  # 波動関数 第2成分 実部
    d = var_list.y[3][i]  # 波動関数 第2成分 虚部
    psi = np.array([[a + b*1j],
                    [c + d*1j]])  # 波動関数

    # 断熱状態に指定
    q_f = eig_vec(var_list.t[i], "plus")  # final state

    PA = np.vdot(q_f, psi)  # probability amplitude
    OP = abs(PA)**2  # ocupation probability
    OP_list.append(OP)

# 終時間における状態0の占有確率
phase = phase_term2 - phase_term1 + phase_term3
P_f_adia = 4 * math.exp(TP) * math.cos(phase)**2
# occupation probability (adiabatic)
P_f_HS = (4 * math.exp(TP) * (1 - math.exp(TP)) * math.cos(phi_s + phase)**2)
# occupation probability (heuristic solution)

# ##出力用プログラム####################################################
# 値を表示する
dic = {
    'v': v,
    'm': m,
    'k': k,
    'P_TLZ': math.exp(TLZ),
    'P_TP': math.exp(TP),
    'P_f_adiabatic': P_f_adia,
    'P_f_HS1': P_f_HS,
    'P_f_num': OP_list[-1],
}
ll = max([len(m) for m in dic.keys()])
for m, v in dic.items():
    print(f'{m:{ll}} : {v}')

# グラフ表示
P_f_adia += t_eval*0
P_f_HS += t_eval*0
plt.plot(t_eval, OP_list)
plt.plot(t_eval, P_f_adia)
plt.plot(t_eval, P_f_HS)
plt.xlim([-3, 3])
plt.ylim([-0.1, 1.1])
plt.show()

# %%
