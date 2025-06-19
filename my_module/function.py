from typing import Callable
import numpy.typing as npt

import numpy as np


def q(time: float, sweep_speed: float) -> float:
    """
    define parameter sweep

    q = adiabatic parameter * time

    Args:
        time (float): time
        sweep_speed (float): sweep speed (adiabatic parameter)

    Returns:
        float: parameter sweep
    """
    return sweep_speed * time


def phi_dot(time: complex, hamiltonian: Callable[..., complex], epsilon: float = 0) -> complex:
    """
    define derivative of azimuthal angle

    Args:
        time (complex): time
        hamiltonian (Callable): Hamiltonian function
        epsilon (float, optional): division by zero regularization parameter. Defaults to 0.

    Returns:
        complex: derivative of azimuthal angle
    """
    H_x: complex = hamiltonian(time, "x")
    H_y: complex = hamiltonian(time, "y")
    H_x_dot: complex = hamiltonian(time, "x_dot")
    H_y_dot: complex = hamiltonian(time, "y_dot")

    num: complex = -H_x * H_y_dot + H_x_dot * H_y
    den: complex = H_x**2 + H_y**2
    return num / (den + epsilon)


def adia_eng(time: complex, Hamiltonian: Callable[..., complex], ut: bool = False, omega: float | None = None) -> complex | float:
    """define adiabatic energy

    Args:
        time (complex): time
        hamiltonian (callable): Hamiltonian
        ut (bool, optional): unitary transformed. Defaults to False.
        F (float, optional): sweep speed. If `ut=True`, this parameter must be specified. Defaults to None.

    Returns:
        complex: adiabatic energy
    """
    H_x: complex = Hamiltonian(time, "x")
    H_y: complex = Hamiltonian(time, "y")
    H_z: complex = Hamiltonian(time, "z")
    phi_d: complex = phi_dot(time, Hamiltonian)

    if ut:
        if omega is None:
            raise ValueError("if `ut` is `True`, the argument 'F' must be specified.")
        else:
            return np.sqrt(H_x**2 + H_y**2 + (H_z + 0.5 * (-omega) * phi_d)**2)
    else:
        return np.sqrt(H_x**2 + H_y**2 + H_z**2)


def eig_vec(time: float, hamiltonian: Callable[..., float], state: str) -> np.ndarray:
    """
    eigenvector

    Args:
        time (float): time
        hamiltonian (Callable): Hamiltonian
        state (str): upper or lower

    Returns:
        array: eigenvector
    """
    H_x: float = hamiltonian(time, "x").real
    H_y: float = hamiltonian(time, "y").real
    H_z: float = hamiltonian(time, "z").real
    adia_energy: float = adia_eng(time, hamiltonian).real

    # lower stateを求めるときは断熱エネルギーを符号反転する
    if state == "lower":
        adia_energy = -adia_energy

    eig_vec: npt.NDArray = np.array([H_x - H_y * 1j, adia_energy - H_z])
    eig_vec /= np.linalg.norm(eig_vec)  # normalization
    return eig_vec


def adia_param(v: float, omega: float, m: float, k: float) -> float:
    return (m - k * v * omega / 4)**2 / (2 * abs(v) * abs(omega))


def TLZ_theoretical(v: float | npt.NDArray, omega: float | npt.NDArray, m: float | npt.NDArray, k: float | npt.NDArray) -> npt.NDArray:
    TLZ = -np.pi * (m + k * v * omega / 4)**2 / (np.abs(v) * np.abs(omega))
    return np.exp(TLZ)


def func_psi_module(time: float, hamiltonian: Callable[..., float], var: list[float], h: float = 1) -> list[float]:
    """
    state vector

    (t_f)における系の状態ベクトル(psi(t_f))を求める関数です。
    psiの第1成分をa+ib, 第2成分をc*idとします。
    var[0]=a,var[1]=b, var[2]=c, var[3]=dとします。

    Args:
        time (float): time
        hamiltonian (Callable): Hamiltonian function
        var (list): 状態ベクトルの各成分を要素とするlist
        h (float, optional): Dirac constant. Default to 1.

    Returns:
        list: 微分方程式
    """

    H_x: float = hamiltonian(time, "x").real
    H_y: float = hamiltonian(time, "y").real
    H_z: float = hamiltonian(time, "z").real

    dadt: float = (+1 / h) * (H_x * var[3] - H_y * var[2] + H_z * var[1])
    dbdt: float = (-1 / h) * (H_x * var[2] + H_y * var[3] + H_z * var[0])
    dcdt: float = (+1 / h) * (H_x * var[1] + H_y * var[0] - H_z * var[3])
    dddt: float = (-1 / h) * (H_x * var[0] - H_y * var[1] - H_z * var[2])

    return [dadt, dbdt, dcdt, dddt]


if __name__ == '__main__':
    print(q(1, 1))
