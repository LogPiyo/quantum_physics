from typing import Callable
import numpy.typing as npt

import numpy as np


def q(t: float, F: float) -> float:
    """
    define parameter sweep

    q = adiabatic parameter * time

    Args:
        t (float): time
        F (float): sweep speed (adiabatic parameter)

    Returns:
        float: parameter sweep
    """
    return -F * t


def phi_dot(t: complex, Ham: Callable[..., complex], eps: float = 0) -> complex:
    """
    define derivative of azimuthal angle

    Args:
        t (complex): time
        eps (float): epsilon

    Returns:
        float: derivative of azimuthal angle
    """
    H_x: complex = Ham(t, "x")
    H_y: complex = Ham(t, "y")
    H_x_dot: complex = Ham(t, "x_dot")
    H_y_dot: complex = Ham(t, "y_dot")

    num: complex = -H_x * H_y_dot + H_x_dot * H_y
    den: complex = H_x**2 + H_y**2
    return num / (den + eps + 1e-10)  # Avoid division by zero with a small epsilon


def adia_eng(t: float, Ham: Callable[..., complex], ut: bool = False, F: float | None = None) -> complex:
    """define adiabatic energy

    Args:
        t (complex): time
        Ham (function): Hamiltonian
        ut (bool, optional): unitary transformed. Defaults to False.
        real (bool, optional): real part. Defaults to False.
        F (float, optional): sweep speed. If `ut=True`, this parameter must be specified. Defaults to None.

    Returns:
        complex: adiabatic energy
    """
    H_x: complex = Ham(t, "x")
    H_y: complex = Ham(t, "y")
    H_z: complex = Ham(t, "z")
    phi_d: complex = phi_dot(t, Ham)

    if ut:
        if F is None:
            raise ValueError("if `ut` is `True`, the argument 'F' must be specified.")
        else:
            return np.sqrt(H_x**2 + H_y**2 + (H_z + 0.5 * (-F) * phi_d)**2)
    else:
        return np.sqrt(H_x**2 + H_y**2 + H_z**2)


def eig_vec(t: float, Ham: Callable[..., float], s: str) -> npt.NDArray:
    """
    eigenvector

    Args:
        t (float): time
        Ham (callable): Hamiltonian
        s (str): upper or lower

    Returns:
        array: eigenvector
    """
    H_x: float = Ham(t, "x").real
    H_y: float = Ham(t, "y").real
    H_z: float = Ham(t, "z").real
    adia_energy: float = adia_eng(t, Ham).real

    # lower stateを求めるときは断熱エネルギーを符号反転する
    if s == "lower":
        adia_energy = -adia_energy

    eig_vec: npt.NDArray = np.array([H_x - H_y * 1j, adia_energy - H_z])
    eig_vec /= np.linalg.norm(eig_vec)  # normalization
    return eig_vec


def adia_param(v: float, F: float, m: float, k: float) -> float:
    return (m - k * v * F / 4)**2 / (2 * abs(v) * abs(F))


def TLZ_theoretical(v: float | npt.NDArray, F: float | npt.NDArray, m: float | npt.NDArray, k: float | npt.NDArray) -> npt.NDArray:
    TLZ = -np.pi * (m + k * v * F / 4)**2 / (np.abs(v) * np.abs(F))
    return np.exp(TLZ)


def func_psi_module(t: float, Ham: Callable[..., float], var: list[float], h: float = 1):
    """
    state vector

    (t_f)における系の状態ベクトル(psi(t_f))を求める関数です。
    psiの第1成分をa+ib, 第2成分をc*idとします。
    var[0]=a,var[1]=b, var[2]=c, var[3]=dとします。

    Args:
        t (float): time
        var (list): 状態ベクトルの各成分を要素とするlist
        h (float, optional): Dirac constant. Default to 1.

    Returns:
        list: 微分方程式
    """

    H_x: float = Ham(t, "x").real
    H_y: float = Ham(t, "y").real
    H_z: float = Ham(t, "z").real

    dadt: float = (+1 / h) * (H_x * var[3] - H_y * var[2] + H_z * var[1])
    dbdt: float = (-1 / h) * (H_x * var[2] + H_y * var[3] + H_z * var[0])
    dcdt: float = (+1 / h) * (H_x * var[1] + H_y * var[0] - H_z * var[3])
    dddt: float = (-1 / h) * (H_x * var[0] - H_y * var[1] - H_z * var[2])

    return [dadt, dbdt, dcdt, dddt]


if __name__ == '__main__':
    print(q(1, 1))
