# 目的
Lim (2015)から導出した遷移確率
```math
P = \exp \left(-4 \frac{(-F)}{|F|} \int_0^{\mathrm{Im} t_j} dv \mathrm{Re} (E|_{t=\mathrm{Re} t_j + iv}) \right)
```
がTakayoshi (2021) 式(14)の遷移確率
```math
\begin{align*}
P(F)
&= \exp \left(-2 \mathrm{Im} \, \int_0^{q_j} \frac{\Delta(q)}{|F|} dq\right) \\
&= \exp \left(-4 \mathrm{Im} \, \int_0^{q_j} \frac{E(q)}{|F|} dq\right)
\end{align*}
```
に一致することを確かめます。
- $`q = -F t`$
- $`E, E(q)`$: 断熱エネルギー
- $`\Delta(q) = 2 E(q)`$
- $`t_j, q_j \in \mathbb{C}`$: 断熱エネルギーの零点

# 方針
Lim (2015) の式(C4)
```math
\int_0^{t_j} dt' (\cdots)
= \int_0^{\mathrm{Re} t_j} du (\cdots) + i \int_0^{\mathrm{Im} t_j} dv (\cdots)|_{t' = \mathrm{Re} t_j + i v}
```
を使います。右辺の第1項は実数，第2項は複素数です。
- $`t' = u + i v \in \mathbb{C}`$: 断熱エネルギーの零点
- $`u, v \in \mathbb{R}`$

この式の両辺の虚部をとると，
```math
\mathrm{Im} \, \int_0^{t_j} E(\cdots) dt
= \int_0^{\mathrm{Im} t_j} dv \mathrm{Re} (\cdots)|_{t=\mathrm{Re} t_j + iv}
```
が成り立ちます。

# 計算
```math
\begin{align*}
    P
    &= \exp \left(-4 \mathrm{Im} \, \int_0^{q_j} \frac{E(q)}{|F|} dq\right) \\
    &= \exp \left(-4 \frac{(-F)}{|F|} \mathrm{Im} \, \int_0^{t_j} E(t) dt\right) \\
    &= \exp \left(-4 \frac{(-F)}{|F|} \int_0^{\mathrm{Im} t_j} dv \mathrm{Re} (E|_{t=\mathrm{Re} t_j + iv}) \right)
\end{align*}
```