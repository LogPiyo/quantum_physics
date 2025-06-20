# self-made libraries
from my_module.function import eig_vec, func_psi_module

# type hint libraries
import numpy.typing as npt
from typing import Callable

# math libraries
import numpy as np
from scipy.integrate import solve_ivp


def calculate_occupation_probability(Hamiltonian: Callable[..., float], initial_time: float, final_time: float, sample_number: int) -> npt.NDArray:
    """
    任意の2準位Hamiltonian系について、
    Shrödinger方程式の数値微分を行うことで、占有確率の時間発展を計算する。
    初期状態がlower stateのとき、upper stateの占有確率を求める。

    Args:
        real_hamiltonian (Callable[..., float]): Hamiltonianを定義する関数
        initial_time (float): 初期時間
        final_time (float): 最終時間

    Returns:
        npt.NDArray: 占有確率の配列
    """
    # solve_ivp計算用の関数
    def func_psi(time, var):
        return func_psi_module(time, Hamiltonian, var)

    # 各Fにおけるpsiの時間発展を計算し，t_fにおけるpsiとFをvar_fに追加する。
    t_eval = np.linspace(initial_time, final_time, sample_number)  # 時間の配列
    initial_state_vecter: list[float] = eig_vec(initial_time, Hamiltonian, "upper").tolist()
    var_init: list[float] = [initial_state_vecter[0].real, initial_state_vecter[0].imag,
                             initial_state_vecter[1].real, initial_state_vecter[1].imag]
    var_list = solve_ivp(func_psi, [initial_time, final_time], var_init, method="LSODA",
                         t_eval=t_eval, rtol=1e-12, atol=1e-12)

    # 占有確率の計算
    state_vecter: npt.NDArray = np.array([[
        var_list.y[0][i] + 1j * var_list.y[1][i],
        var_list.y[2][i] + 1j * var_list.y[3][i]
    ] for i in range(sample_number)])
    final_state_vecter: npt.NDArray = np.array([eig_vec(var_list.t[i], Hamiltonian, "lower") for i in range(sample_number)])
    transition_amplitude: npt.NDArray = np.sum(np.conj(final_state_vecter) * state_vecter, axis=1)
    occupation_probability_arr: npt.NDArray = np.abs(transition_amplitude)**2

    return occupation_probability_arr


if __name__ == "__main__":
    # Example usage
    def example_hamiltonian(t, component):
        # Define a simple Hamiltonian for demonstration
        return t if component == 'x' else 0

    initial_time = 0.0
    final_time = 10.0
    sample_number = 100

    occupation_probabilities = calculate_occupation_probability(example_hamiltonian, initial_time, final_time, sample_number)
    print(occupation_probabilities)
