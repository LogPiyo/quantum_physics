## Zenerの証明法
Lanau-Zener公式をZenerの方法で証明します。

```math
H = 
\begin{pmatrix}
    -\frac{vt}{2} & \Delta_0 \\
    \Delta_0 & \frac{vt}{2}
\end{pmatrix}
```
とします。また，系の状態ベクトルは，
```math
  |\psi(t) \rangle = C_1(t) |1\rangle + C_2(t) |2\rangle
```
です。さらに，系の初期状態は，$`C_1(-\infty) = 1, C_2(-\infty) = 0`$としておきます。このとき，Shrödinger方程式より，
```math
\begin{align}
    i\frac{d}{dt} C_1(t)
    &= -\frac{vt}{2} C_1(t) + \Delta_0 C_2(t) \\
    i\frac{d}{dt} C_2(t)
    &= \Delta_0 C_1(t) + \frac{vt}{2} C_2(t) \\
\end{align}
```
となります。式(3)を微分して式(4)を代入すると，
```math
\begin{align}
    i\frac{d^2}{dt^2} C_1
    &= -\frac{v}{2} C_1(t) -\frac{vt}{2} C_1^{\prime} + \Delta_0 C_2^{\prime}\\
    &= -\frac{v}{2} C_1(t) -\frac{vt}{2} C_1^{\prime} + \Delta_0 \left( \frac{\Delta_0}{i} C_1 + \frac{vt}{2i} C_2 \right)
\end{align}
```
です。また，式(4)を微分して式(3)を代入すると。同様の計算から，
```math
\begin{align}
    i\frac{d^2}{dt^2} C_2
    = \frac{v}{2} C_2 + \frac{vt}{2} C_2^{\prime} + \Delta_0 \left( \frac{\Delta_0}{i} C_2 - \frac{vt}{2i} C_1 \right)
\end{align}
```
が得られます。したがって，式(7)および式(8)を整理すると，
```math
\begin{align}
    \frac{d^2 C_1}{dt^2} + \left( \Delta_0^2 - \frac{iv}{2} \right) C_1 + \frac{v^2 t^2}{4} C_1
    = 0 \\
    \frac{d^2 C_2}{dt^2} + \left( \Delta_0^2 + \frac{iv}{2} \right) C_2 + \frac{v^2 t^2}{4} C_2
    = 0 \\
\end{align}
```
という$`C_1`$と$`C_2`$について独立な2つの方程式が得られます。ここで，
```math
\begin{align}
    z &= i\sqrt{v} e^{i\frac{\pi}{4}} t \\
    n &= i\delta \\
    \delta &= \frac{\Delta_0^2}{v} \\
    C_1 &= w_1(z) \\
    C_2 &= w_2(z)
\end{align}
```
という変数変換を行うと，
```math
\begin{align}
    \frac{d^2 w_1}{d z^2} + \left(n + \frac{1}{2} - \frac{z^2}{4}\right) w_1
    = 0 \\
    \frac{d^2 C_2}{d z^2} + \left(n - \frac{1}{2} - \frac{z^2}{4}\right) w_2
    = 0
\end{align}
```
が得られます。これらはWeberの微分方程式です。ここで，Weberの微分方程式の基本解である放物柱関数$`D_n(z)`$を用いて，
```math
\begin{align}
    w_1 &\propto D_n(z)\\
    w_2 &\propto D_{n-1}(z)
\end{align}
```
として構いません。なぜなら，式(\ref{z_def})より，$`t<0`$で，
```math
z
= i\sqrt{v} e^{i\frac{\pi}{4}} |t|
= \sqrt{v} e^{-i\frac{\pi}{4}} |t|
```
と書けるため，
```math
|\arg z| = \frac{\pi}{4} < \frac{3\pi}{4}
```
であり，放物柱関数$`D_n(z)`$に成り立つ定理
```math
|z| \gg 1のとき，|\arg z| < \frac{3\pi}{4}ならば
D_n(z) \rightarrow e^{-\frac{z^2}{4}} z^n
```
を用いることができるからです。したがって，$`t \rightarrow - \infty`$で，
```math
\begin{align}
    D_n(z) &= e^{-\frac{z^2}{4}} z^n\\
    D_{n-1}(z) &= e^{-\frac{z^2}{4}} z^{n-1}
\end{align}
```
と書けます。ここで，
```math
\begin{align}
    z^n
    &= z^{i\delta}\\
    &= \exp(i\delta \log z)\\
    &= \exp(i\delta (\log \sqrt{v} e^{-\frac{i\pi}{4}} |t|))\\
    &= \exp \left( i\delta (\log (\sqrt{v} |t|) - \frac{i\pi}{4} \right)
\end{align}
```
より，$`x \in \mathbb{R}`$のとき，
```math
\begin{align}
    D_n(z)
    = \exp \left( i\delta (\log (\sqrt{v} |t|) + \frac{\pi \delta}{4} \right) \\
    D_{n-1}(z)
    = \frac{\exp \left( i\delta (\log (\sqrt{v} |t|) + \frac{\pi \delta}{4} \right)}{\sqrt{v} e^{-\frac{i\pi}{4}} |t|}
\end{align}
```
であるため，$`t \rightarrow -\infty`$で，
```math
\begin{align}
    |D_n(z)| &\rightarrow \text{有限} \\
    |D_{n-1}(z)| &\rightarrow 0
\end{align}
```
となります。これは，系の初期条件
```math
\begin{align}
    C_1(-\infty) \ne 0\\
    C_2(-\infty) = 0
\end{align}
```
を満たします。したがって，任意定数$`A,B`$を用いて，
```math
\begin{align}
    C_1(t) &= w_1(z) = A D_n(z),\\
    C_2(t) &= w_2(z) = B D_{n-1}(z)
\end{align}
```
として構いません。したがって，式(\ref{D_n})および式(\ref{D_n-1})より，$`t \rightarrow -\infty`$のとき，
```math
\begin{align}
    \lim_{t \rightarrow -\infty} C_1(t)
    &= A \exp \left( \frac{ivt^2}{4} + \frac{\pi \delta}{4} + i\delta \log \sqrt{v} |t| \right) \\
    \lim_{t \rightarrow -\infty} C_2(t)
    &= 0
\end{align}
```
です。
一方，$`t \rightarrow \infty`$のとき，
```math
\begin{align}
    z &= \sqrt{v} e^{i\frac{3\pi}{4}},\\
    z^n
    &= z^{i\delta}\\
    &= \exp(i\delta \log z)\\
    &= \exp ( i\delta \log (\sqrt{v} t e^{i\frac{3\pi}{4}}))\\
    &= \exp \left( i\delta (\log \sqrt{v} t) + i\frac{3\pi}{4} \right)
\end{align}
```
となります。ここで，放物柱関数$`D_n(z)`$に成り立つ定理
```math
\frac{\pi}{4} < \arg z < \frac{3\pi}{4}のとき，
D_n(z) \rightarrow e^{-\frac{z^2}{4}} z^n - \frac{\sqrt{2\pi}}{\Gamma(-n)} e^{in\pi} e^{\frac{z^2}{4}} z^{-n-1}
```
を用います。このとき，$`C_1(t)`$については，$`t \rightarrow \infty`$で(第2項)$`\rightarrow 0`$より，
```math
\lim_{t\rightarrow \infty} C_1(t)
= A \exp \left( \frac{ivt^2}{4} - \frac{e\pi \delta}{4} + i\delta \log (\sqrt{v}t) \right)
```
となります。また，$`C_2(t)`$については，$`t \rightarrow \infty`$(第1項)$`\rightarrow 0`$より，
```math
\begin{align}
    \lim_{t \rightarrow \infty} C_2(t)
    &= -B \frac{\sqrt{2 \pi}}{\Gamma(n+1)} e^{i(n-1) \pi} e^{\frac{z^2}{4}} z^{-n} \\
    &= -B \frac{\sqrt{2 \pi}}{\Gamma(1-i \delta)} \exp \left(i(i \delta-1) \pi-\frac{i v t^2}{4}-i \delta \log (\sqrt{v} t) \frac{3 \pi   \delta}{4}\right) \\
    & = B \frac{\sqrt{2 \pi}}{\Gamma(1-i \delta)} \exp \left(-\frac{\pi \delta}{4}-\frac{i v t^2}{4}-i \delta \log \sqrt{v} t\right)
\end{align}
```
となります。したがって，式(\ref{C1_inf})および式(\ref{C1_lim2})より，
```math
\frac{\lim_{t \rightarrow \infty} C_1(t)}{\lim_{t \rightarrow -\infty} C_1(t)}
= \exp (-\pi \delta)
```
であり，Landau-Zener公式
```math
P
= \left| \frac{\lim_{t \rightarrow \infty} C_1(t)}{\lim_{t \rightarrow -\infty} C_1(t)} \right|^2
= \exp(-2\pi \delta)
```
が示されました。