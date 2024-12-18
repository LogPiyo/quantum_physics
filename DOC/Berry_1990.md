# Berry (1990) のレビュー
## 概要

$$
\begin{align*}
    H(t)
    &=
    \begin{pmatrix}
        B_z & B_x - i B_y \\
        B_x + i B_y & B_z
    \end{pmatrix} \\
    &= H
    \begin{pmatrix}
        \cos{\theta} & \sin{\theta} e^{-i \phi} \\
        \sin{\theta} e^{-i \phi} & \cos{\theta}
    \end{pmatrix} \\
\end{align*}
$$
というHamlitonianにおいて，遷移確率は
$$
\begin{align*}
    P
    &= |\langle u_-(+\infty) | \psi(+\infty) \rangle|^2 \\
    &\approx \exp{(-\Gamma_d)} \exp{(\Gamma_g)}
\end{align*}
$$
で与えられる。ただし，
$$
\begin{align*}
    \Gamma_d
    &= \frac{4}{\hbar \delta} \mathrm{Im} \int_0^{\tau_c} d\tau |H(\tau)| \\
    \Gamma_g
    &= 2 \mathrm{Im} \int_0^{\tau_c} d\tau \dot{\phi} \cos{\theta} \\
    &= -\frac{1}{2} \mathrm{Im} \oint d\tau \dot{\phi} \cos{\theta}
\end{align*}
$$
である。

