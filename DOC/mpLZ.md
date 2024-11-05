# multiple-passage Landau-Zenerモデル

Hamiltonian
```math
H_\mathrm{MTLZ}(t)
=
\begin{pmatrix}
    \Delta_z \sin \omega t & -\varepsilon_0 \cos \omega t\\
    -\varepsilon_0 \cos \omega t & -\Delta_z \sin \omega t
\end{pmatrix}
```
をmultiple-passage Landau-Zenerモデルと呼びます。$`n`$回交差してlower-stateへ遷移する確率は，
```math
\left\{
\begin{align*}
    P_{2m-1}
    &= 1-q \left( \frac{\sin (2m-1) \xi}{\sin \xi} \right)^2 \\
    P_{2m} 
    &= q \left( \frac{\sin (2m) \xi}{\sin \xi} \right)^2
\end{align*}
\right.
```
で与えられます。

## 遷移確率の導出法

$C(t) =\sin \omega t$のとき，Hamiltonianは，
```math
\hat{H}
=
\begin{pmatrix}
    \varepsilon_0 \cos \omega t & \Delta_0 \sin \omega t\\
    \Delta_0 \sin \omega t & -\varepsilon_0 \cos \omega t
\end{pmatrix}
```
である．また，$M$および$G$に対応する行列は，それぞれの複素共役行列である．そのため，初期状態を$\psi_0 = (1,0)$とすると，
```math
\begin{align}
    \psi_1 = G^{\frac{1}{2}} M_1 G^{*\frac{1}{2}} \psi_0\\
    \psi_2 = G^{*\frac{1}{2}} M_2 G^{\frac{1}{2}} \psi_0
\end{align}
```
である．ここで，表式を簡単にするため，
```math
\begin{align}
    T &= G^{\frac{1}{2}} M G^{*\frac{1}{2}}\\
    T^* &= G^{*\frac{1}{2}} M^* G^{\frac{1}{2}}
\end{align}
```
としよう．このとき，
```math
\begin{align}
    \psi_{2m-1}
    &= T (T^* T)^{m-1} \psi_0 \\
    \psi_{2m}
    &= (T^* T)^{m} \psi_0
\end{align}
```
と書ける．また，
```math
\begin{align}
    S
    = 
    \begin{pmatrix}
        \sqrt{1-q} e^{-i(\phi_s + \theta)} & \sqrt{q}\\
        -\sqrt{q} & \sqrt{1-q} e^{i(\phi_s + \theta)}
    \end{pmatrix}
\end{align}
```
を定義すると，
```math
T^* T = - S^2
```
という関係式が成り立つことが確かめられる．さらに，
```math
\begin{align}
    S^2
    &=
    \begin{pmatrix}
        \sqrt{1-q} e^{-i(\phi_s + \theta)} & \sqrt{q}\\
        -\sqrt{q} & \sqrt{1-q} e^{-i(\phi_s + \theta)}
    \end{pmatrix}
    \begin{pmatrix}
        \sqrt{1-q} e^{-i(\phi_s + \theta)} & \sqrt{q}\\
        -\sqrt{q} & \sqrt{1-q} e^{-i(\phi_s + \theta)}
    \end{pmatrix} \\
    &=
    \begin{pmatrix}
        \sqrt{1-q} e^{-i(\phi_s + \theta)} \cos \xi -1 & 2 \sqrt{q} \cos \xi\\
        -2 \sqrt{q} \cos \xi & \sqrt{1-q} e^{i(\phi_s + \theta)} \cos \xi -1
    \end{pmatrix} \\
    &=
    2 \cos \xi S - 1
\end{align}
```
が成り立つ．ただし，
```math
\cos \xi = \sqrt{1-q} \cos (\phi_s + \theta)
```
である．特に，$\cos \xi = 0$のとき，$\psi_{2m} = \psi_0$となるから，完全に初期状態に戻る．このときの条件は，
```math
\phi_s + \theta = \left(k + \frac{1}{2}\right) \pi \quad (k = 0,1,2,\ldots)
```
で与えられる．ここで，
```math
S^n = \alpha_n + \beta_n S
```
を仮定すると，
```math
\begin{align}
    &S^2 - 2 \cos\xi S+ 1 = 0\\
    \Leftrightarrow \quad &S^{n+1} - 2 \cos \xi S^n + S^{n-1}= 0\\
    \Leftrightarrow \quad &(\alpha_{n+1} + \beta_{n+1} S) - 2\cos \xi (\alpha_n + \beta_n S) + (\alpha_{n-1} + \beta_{n-1} S) = 0
\end{align}
```
より，$\alpha_{n+1} = - \beta_n$とすると，
```math
\beta_{n+1} - 2 \cos\xi \beta_n + \beta_{n-1} = 0
```
という2階差分方程式を得る．式(\ref{DE})を解くためには，$\beta_n = \rho^n$と置けばよい．すると，式(\ref{DE})は，
```math
\begin{equation}
  \rho^{n+1} - 2 \cos \xi \rho^n + \rho^{n-1} = 0
\end{equation}
```
となるから，
```math
\begin{equation}
  \rho = \cos \xi \pm i \sin \xi = e^{\pm i \xi}
\end{equation}
```
と決まる．したがって，$B_n$の一般解は，任意定数$C_1$，$C_2$を用いて，
```math
\begin{equation}
  \beta_n = C_1 e^{in\xi} + C_2 e^{-in\xi}
\end{equation}
```
と書ける．式(\ref{relation_S})から明らかな条件$\beta_0 = 0$，$\beta_1=0$を用いると，
```math
\begin{equation}
  \left\{ \,
  \begin{aligned}
    & C_1 + C_2 &= 0 \\
    & C_1 e^{i\xi} + C_2 e^{-i\xi} &= 1 
  \end{aligned}
  \right.
\end{equation}
```
という2元連立方程式が得られる．よって，
```math
\begin{equation}
  \left\{ \,
    \begin{aligned}
    &  C_1 &= \frac{1}{2i\sin \xi}\\
    &  C_2 &= -\frac{1}{2i\sin \xi}
    \end{aligned}
  \right.
\end{equation}
```
と決定される．これらを，式(\ref{DE})に代入すると，
```math
\begin{align}
  \beta_n
  &= \frac{1}{2i\sin \xi} (e^{in\xi} - e^{-in\xi})\\
  &= \frac{\sin n \xi}{\sin \xi}
\end{align}
```
となる．したがって，式(\ref{psi_2m-1})は，
```math
\begin{align}
  \psi_{2m-1}
  &= T (T^* T)^{m-1} \psi_0\\
  &= T (-1)^{m-1} S^{2m-2} \psi_0\\
  &= 
  \begin{pmatrix}
    \sqrt{q} \beta_{2m-1}\\
    \sqrt{1-q} e^{-i(\phi_s+\theta)} \beta_{2m-1} - \beta_{2m-2}
  \end{pmatrix}
\end{align}
```
となる．また，同様の計算から，式(\ref{psi_2m})は，
```math
\begin{align}
  \psi_{2m} = 
  \begin{pmatrix}
    \sqrt{1-q} e^{-i(\phi_s+\theta)} \beta_{2m} - \beta_{2m-1}\\
    -\sqrt{q} \beta_{2m}
  \end{pmatrix}
\end{align}
```
となる．したがって，この系が，$t = \frac{(2m-1)\pi}{\omega}, \frac{2m \pi}{\omega}$においてそれぞれ状態$|1\rangle$にある確率$P_{2m-1}, P_{2m}$は，
```math
\begin{align}
  P_{2m-1}
  &= 1-q |\beta_{2m-1}|^2\\
  &= 1-q \left( \frac{\sin (2m-1) \xi}{\sin \xi} \right)^2,\\
  P_{2m} 
  &= q |\beta_{2m}|^2\\
  &= q \left( \frac{\sin (2m) \xi}{\sin \xi} \right)^2
\end{align}
```
である．