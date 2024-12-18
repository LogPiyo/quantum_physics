# Majoranaの証明法
Landau-Zener公式をMajoranaの方法で証明します。
## Hamiltonian
```math
H
=
\begin{pmatrix}
-\frac{1}{2} \nu t & -\frac{1}{2} \Delta \\
-\frac{1}{2} \Delta & \frac{1}{2} \nu t\\
\end{pmatrix}
```

## 計算
Shrödinger方程式より，
```math
i \hbar \frac{\partial}{\partial t} | \psi \rangle
= -\frac{1}{2} (\Delta \sigma_x + \nu t \sigma_z) | \psi \rangle
```
である。ここで，
```math
\tau = \sqrt{\frac{\nu}{2 \hbar}} t \, , \delta = \frac{\Delta^2}{4 \nu \hbar}
```
というパラメータを導入する。また，
```math
\alpha = f \exp{\left(\frac{1}{2} \tau^2 \right)}\, , \,
\beta = g \exp{\left(-\frac{i}{2} \tau^2 \right)}
```
とする。このとき，
```math
\left\{
    \begin{align*}
        i \hbar \frac{\partial}{\partial t} \alpha
        &= -\frac{1}{2} (\nu t \alpha + \Delta \beta) \\
        i \hbar \frac{\partial}{\partial t} \beta
        &= -\frac{1}{2} (\Delta \alpha - \nu t \beta)
    \end{align*}
\right.
```
$`\Leftrightarrow`$
```math
\left\{
    \begin{align*}
        i \hbar \sqrt{\frac{\nu}{2 \hbar}} \dot{\alpha}
        &= -\frac{1}{2} (\nu t \alpha + \Delta \beta) \\
        i \hbar \sqrt{\frac{\nu}{2 \hbar}} \dot{\beta}
        &= -\frac{1}{2} (\Delta \alpha - \nu t \beta)
    \end{align*}
\right.
```
$`\Leftrightarrow`$
```math
\left\{
    \begin{align*}
        i \hbar \sqrt{\frac{\nu}{2 \hbar}} \dot{\alpha}
        &= -\frac{1}{2} (\nu \sqrt{\frac{2 \hbar}{\nu}} \tau \alpha + \Delta \beta) \\
        i \hbar \sqrt{\frac{\nu}{2 \hbar}} \dot{\beta}
        &= -\frac{1}{2} (\Delta \alpha - \sqrt{\frac{2 \hbar}{\nu}} \tau \beta)
    \end{align*}
\right.
```
$`\Leftrightarrow`$
```math
\left\{
    \begin{align*}
        i \hbar \sqrt{\frac{\hbar \nu}{2}} \dot{\alpha}
        &= -\frac{1}{2} (\sqrt{2 \hbar \nu} \tau \alpha + \Delta \beta) \\
        i \hbar \sqrt{\frac{\hbar \nu}{2}} \dot{\beta}
        &= -\frac{1}{2} (\Delta \alpha - \sqrt{2 \hbar \nu} \tau \beta)
    \end{align*}
\right.
```
$`\Leftrightarrow`$
```math
\left\{
    \begin{align*}
        i \hbar \sqrt{\frac{1}{2}} \dot{\alpha}
        &= -\frac{1}{2} (\sqrt{2} \tau \alpha + \frac{\Delta}{\sqrt{\hbar \nu}} \beta) \\
        i \hbar \sqrt{\frac{1}{2}} \dot{\beta}
        &= -\frac{1}{2} (\frac{\Delta}{\sqrt{\hbar \nu}} \alpha - \sqrt{2} \tau \beta)
    \end{align*}
\right.
```
$`\Leftrightarrow`$
```math
\left\{
    \begin{align*}
        i \hbar \sqrt{\frac{1}{2}} \dot{\alpha}
        &= -\frac{1}{2} (\sqrt{2} \tau \alpha + 2 \sqrt{\delta}  \beta) \\
        i \hbar \sqrt{\frac{1}{2}} \dot{\beta}
        &= -\frac{1}{2} (2 \sqrt{\delta} \alpha - \sqrt{2} \tau \beta)
    \end{align*}
\right.
```
$`\Leftrightarrow`$
```math
\left\{
    \begin{align*}
        i \dot{\alpha}
        &= -(\tau \alpha + \sqrt{2 \delta}  \beta) \\
        i \dot{\beta}
        &= -(\sqrt{2 \delta} \alpha - \tau \beta)
    \end{align*}
\right.
```
です。したがって，
```math
\left\{
    \begin{align*}
        i (\dot{f} \exp(\frac{i}{2} \tau^2) + i \tau f \exp(\frac{i}{2} \tau^2))
        &= -\tau f \exp(\frac{i}{2} \tau^2) - \sqrt{2 \delta} g \exp(-\frac{i}{2} \tau^2) \\
        i (\dot{g} \exp(-\frac{i}{2} \tau^2) - i \tau g \exp(-\frac{i}{2} \tau^2))
        &= \sqrt{2 \delta} f \exp(\frac{i}{2} \tau^2) + \tau g \exp(-\frac{i}{2} \tau^2)
    \end{align*}
\right.
```
$`\Leftrightarrow`$
```math
\left\{
    \begin{align*}
        \dot{f}
        &= i \sqrt{2 \delta} g \exp(-i \tau^2) \\
        \dot{g}
        &= i \sqrt{2 \delta} f \exp(i \tau^2)
    \end{align*}
\right.
```
です。すなわち，
```math
\begin{align*}
    \ddot{f}
    &= i \sqrt{2 \delta} (\dot{g} \exp(-i \tau^2) + g (-2 i \tau) \exp(-i \tau^2)) \\
    &= i \sqrt{2 \delta} (i \sqrt{2 \delta} f + \frac{1}{i \sqrt{2 \delta}} (-2 i \tau) \dot{f}) \\
    &= -2 \delta f - 2 i \tau \dot{f}
\end{align*}
```
より，
```math
\ddot{f} + 2 i \tau \dot{f} + 2 \delta f = 0
```
が成り立ちます。$`g`$についても同様に，
```math
\ddot{g} - 2 i \tau \dot{g} + 2 \delta g = 0
```
が成り立ちます。この方程式を解くために，両側Laplace変換
```math
\mathcal{L}[f(\tau)]
= \int_{-\infty}^{\infty} e^{-s \tau} f(\tau) d\tau
= F(s)
```
を導入します。このとき，
```math
\mathcal{L}[\ddot{f}] + 2 i \mathcal{\tau \dot{f}} + 2 \delta \mathcal{l}[f]
= 0
```
と変形できます。ここで，
```math
\begin{align*}
    \mathcal{L}[\ddot{f}]
    &= s^2 F(s) \\
    \mathcal{L}[\tau \dot{f}]
    &= -F'(s) \\
    \mathcal{L}[\tau \dot{f}]
    &= \int_{\infty}^{\infty} e^{-s \tau} \tau \dot{f} d\tau \\
    &= e^{-s \tau} \tau f|_{-\infty}^{\infty} - \int_{\infty}^{\infty} (-s) e^{-s \tau} \tau f d\tau - \int_{\infty}^{\infty} e^{-s \tau} f d\tau \\
    &= s \mathcal{L}[\tau f] - \mathcal{L}[f] \\
    &= -(s F'(s) + F(s))
\end{align*}
```
より，
```math
s^2 F(s) - 2 i (F(s) + s F'(s)) + 2 \delta F(s) = 0
```
を得ます。この解は，
```math
F(s) = C_{\delta} \exp(-\frac{i s^2}{4}) s^{-1 - i \delta}
```
です (要追記)。ここで，Laplace逆変換を行います。
そのためには，鞍点法を使います。

## 鞍点法
```math
I(\tau)
= \int_L \phi(z) e^{\tau^2 \mu(z)} dz
```
という積分に対して，3つの条件

1. $`\mathrm{Re} = \max \mathrm{Re}[\mu(z)]`$となる$`z_0`$が存在する
1. $`\int_L |\phi(z)| \exp(\tau_0 \mathrm{Re} \mu(z)) dz < \infty`$を満たす
1. $`\mu'(z_0) = 0, \mu''(z_0) \ne 0, \frac{d^2}{dy^2} (\mathrm{Re} \mu(z_0 + y \lambda)) = \mathrm{const.} < 0`$を満たす

が成立するとき，
```math
I(\tau) \approx \sqrt{\frac{2 \pi}{\tau^2}} \frac{\phi(z) e^{\tau^2 \mu(z_0)}}{\sqrt{-\mu''(z_0)}}
```
となることが知られています。

## 計算続き

今，$`s = \tau z, ds = \tau dz`$とおくと，
```math
f(\tau)
= C_{\delta} \tau^{-i \delta} \int_L \exp(\tau^2 (z - i \frac{z^2}{4})) z^{-(i \delta + 1)} dz
```
であるから (要追記)，
```math
f(\tau) = C_{\delta} \tau^{-i \delta} I(\tau)
```
として比較すると，
```math
\left\{
    \begin{align*}
        \phi(z) 
        = z^{-(i \delta + 1)} \\
        \mu(z)
        = z - i \frac{z^2}{4}
    \end{align*}
\right.
```
です。条件3から，
```math
\mu'(z) = 1 - \frac{i z_0}{2} = 0 \quad \therefore z_0 = -2 i
```
より，鞍点付近で，
```math
\begin{align*}
    f(\tau)
    &= C_{\delta} \tau^{-i \delta} \sqrt{\frac{2 \pi}{\tau^2}} \frac{(-2 i)^{-(i \delta + 1) e^{\tau^2 (-2 i + \frac{i}{4} \cdot 4)}}}{\sqrt{\frac{i}{2}}} \\
    &= C_{\delta} \sqrt{4 \pi} (-2 i \tau)^{-1 - i \delta} \exp(-i \tau^2 + i \frac{3}{4} \pi)
\end{align*}
```
です。また，ゼロ付近では，$`s^2`$の項を無視して，$`x = s \tau`$とすると，
```math
\begin{align*}
f(\tau)
&= C_{\delta} \int_{L_0} e^x \tau^{1 + i \delta} x^{-1 - i \delta} \tau^{-1} dx \\
&= C_{\delta} \tau^{i \delta} \int_{L_0} e^x x^{-1 - i \delta} dx \\
&\approx C_{\delta} \tau^{i \delta} \frac{2 \pi i}{\Gamma(1 + i \delta)}
\end{align*}
```
となります。したがって，$`\tau < 0`$で，
```math
\left\{
    \begin{align*}
        \alpha(\tau)
        &= C_{\delta} \sqrt{4 \pi} (-2 i \tau)^{-1 - i \delta} \exp(-i \frac{\tau^2}{2} + i \frac{3}{4} \pi) \\
        \beta(\tau)
        &= C_{\delta} \sqrt{\frac{2 \pi}{\delta}} (-2 i \tau)^{-i \delta} \exp(-i \tau^2 + i \frac{3}{4} \pi)
    \end{align*}
\right.
```
$`\tau > 0`$で，
```math
\left\{
    \begin{align*}
        \alpha(\tau)
        &= C_{\delta} \sqrt{4 \pi} (-2 i \tau)^{-1 - i \delta} \exp(-i \frac{\tau^2}{2} + i \frac{3}{4} \pi) + C_{\delta} \tau^{i \delta} \frac{2 \pi i}{\Gamma(1 + i \delta)} \exp(\frac{i \tau^2}{2}) \\
        \beta(\tau)
        &= C_{\delta} \sqrt{\frac{2 \pi}{\delta}} (-2 i \tau)^{-i \delta} \exp(-i \tau^2 + i \frac{3}{4} \pi) + C_{\delta} \sqrt{\frac{\delta}{2}} \tau^{i \delta - 1} \frac{2 \pi i}{\Gamma(1 + i \delta)} \exp(\frac{i \tau^2}{2}) \, ,
    \end{align*}
\right.
```
が成り立ちます。$`|\alpha|^2 = 0, |\beta|^2 = 0`$とすると，
```math
C_{\delta} = \sqrt{\frac{\delta}{2 \pi}} \exp(-\frac{\pi \delta}{2})
```
と決まります。ここで，$`\tau^{-1}`$の項を無視すると，$`\tau \to -\infty`$のとき，
```math
\left\{
    \begin{align*}
        \alpha
        &\to 0 \\
        \beta
        &\to (-2 i \tau)^{-i \delta} \exp(\frac{i \pi}{4} - \frac{\pi \delta}{2} - \frac{i \tau^2}{2}) \, ,
    \end{align*}
\right.
```
$`\tau \to \infty`$のとき，
```math
\left\{
    \begin{align*}
        \alpha
        &\to \tau^{i \delta} \frac{\sqrt{2 \pi \delta}}{\Gamma(1 + i \delta)} \exp(\frac{i \pi}{2} - \frac{\pi \delta}{2} - \frac{i \tau^2}{2}) \\
        \beta
        &\to (-2 i \tau)^{-i \delta} \exp(\frac{i \pi}{4} - \frac{\pi \delta}{2} - \frac{i \tau^2}{2}) \, ,
    \end{align*}
\right.
```
です。ここで，$`\tau > 0`$のとき，
```math
\beta(\tau \to \infty)
= (2 \tau)^{-i \delta} \exp(-\pi \delta + \frac{i \pi}{4} - \frac{i \tau^2}{2})
```
より，
```math
\mathcal{P}
= |\beta(\tau \to \infty)|^2
= \exp(-2 \pi \delta)
```
となり，Landau-Zener公式が示されました。