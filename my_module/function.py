import numpy as np


def q(t, F):
    """
    define parameter sweep

    q = adiabatic parameter * time

    Args:
        t (float): time
        F (float): sweep speed

    Returns:
        float: parametersweep
    """
    return -F * t


def phi_dot(t, Ham, eps=0):
    """
    define derivative of azimuthal angle

    Args:
        t (complex): time

    Returns:
        float: derivative of azimuthal angle
    """
    num = -Ham(t, "x") * Ham(t, "y_dot") + Ham(t, "x_dot") * Ham(t, "y")
    den = Ham(t, "x")**2 + Ham(t, "y")**2
    return num / (den + eps)


def adia_eng(t, Ham, ut=False, real=False, F=None):
    """define adiabatic energy

    Args:
        t (complex): time
        Ham (function): Hamiltonian
        ut (bool, optional): unitary transformed. Defaults to False.
        real (bool, optional): real part. Defaults to False.

    Returns:
        complex: adiabatic energy
    """
    Ham_x = Ham(t, "x")
    Ham_y = Ham(t, "y")
    Ham_z = Ham(t, "z")
    phi_d = phi_dot(t, Ham)

    if ut:
        return np.sqrt(Ham_x**2 + Ham_y**2 + (Ham_z + 0.5*(-F)*phi_d)**2)
    elif real:
        return np.sqrt(Ham_x**2 + Ham_y**2 + Ham_z**2)
    return np.sqrt(Ham_x**2 + Ham_y**2 + Ham_z**2)


def eig_vec(t, Ham, s):
    """
    eigenvector

    Args:
        t (float): time
        Ham (function): Hamiltonian
        s (state): upper or lower

    Returns:
        array: eigenvector
    """
    energy = adia_eng(t, Ham, real=True)

    # lower stateを求めるときは断熱エネルギーを符号反転する
    if s == "lower":
        energy = -energy

    eig_vec = np.array([Ham(t, "x", real=True) - Ham(t, "y", real=True) * 1j, energy - Ham(t, "z", real=True)])
    eig_vec /= np.linalg.norm(eig_vec)  # normalization
    return eig_vec


def adia_param(v, F, m, k):
    return (m - k * v * F / 4)**2 / (2 * abs(v) * abs(F))


def TLZ_theoretical(v, F, m, k):
    TLZ = -np.pi * (m + k*v*F/4)**2 / (abs(v) * abs(F))
    return np.exp(TLZ)


def func_psi_module(t, Ham, var, h=1):
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
    dadt = (1/h) * (Ham(t, "x")*var[3] - Ham(t, "y")*var[2] + Ham(t, "z")*var[1])
    dbdt = (-1/h) * (Ham(t, "x")*var[2] + Ham(t, "y")*var[3] + Ham(t, "z")*var[0])
    dcdt = (1/h) * (Ham(t, "x")*var[1] + Ham(t, "y")*var[0] - Ham(t, "z")*var[3])
    dddt = (-1/h) * (Ham(t, "x")*var[0] - Ham(t, "y")*var[1] - Ham(t, "z")*var[2])

    return [dadt, dbdt, dcdt, dddt]


if __name__ == '__main__':
    print(q(1, 1))
