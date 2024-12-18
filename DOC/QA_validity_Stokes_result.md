# Stokes位相のグラフの妥当性
## 前提
### 断熱パラメータ
Landau-Zenerモデルには，
$$
\begin{align*}
    \delta_{\mathrm{LZ}}
    &= \frac{\Delta_z^2}{2 |\varepsilon_0| |\omega|} \, ,
\end{align*}
$$
twisted Landau-Zenerモデルには，
$$
\begin{align*}
    \bar{\delta}_{\mathrm{TLZ}}
    &\approx \frac{(\Delta_z +\frac{\Delta_y \omega}{\varepsilon_0})^2}{2 |\varepsilon_0| |\omega|}
\end{align*}
$$
を用いる。

#### 近似の意味
上の断熱パラメータの表式は，線形なエネルギー変化を前提としている。ユニタリ変換後のtwisted Landau-Zenerモデル
$$
\begin{pmatrix}
    m + \frac{\kappa_{||} \nu F}{4} & \sqrt{(\nu F t)^2 + (\frac{1}{2} \kappa_{||} \nu^2 (F t)^2)^2}\\
    \sqrt{(\nu F t)^2 + (\frac{1}{2} \kappa_{||} \nu^2 (F t)^2)^2} & - (m + \frac{\kappa_{||} \nu F}{4})
\end{pmatrix}
$$
が線形なエネルギー変化と見なせるためには，
$$
\begin{align*}
    \nu F t &\gg \frac{1}{2} \kappa_{||} \nu^2 (F t)^2 \\
    \Leftrightarrow 
    1 &\gg \frac{1}{2} \kappa_{||} \nu F t    
\end{align*}
$$
を満たす必要がある。このとき，
$$
\begin{align*}
    \sqrt{(\nu F t)^2 + \left(\frac{1}{2} \kappa_{||} \nu^2 (F t)^2 \right)^2}
    &= \sqrt{(\nu F t)^2 \left(1 + \left(\frac{1}{2} \kappa_{||} \nu F t \right)^2 \right)} \\
    &\approx |\nu F t| \left(1 + \frac{1}{2} \left(\frac{1}{2} \kappa_{||} \nu F t \right)^2 \right) \\
\end{align*}
$$
である。    
### Stokes位相
$$
\begin{align*}
    \varphi_s(\delta)
    = \frac{\pi}{4} + \delta (\ln \delta -1) + \mathrm{Arg}\, \Gamma (1 - i\delta)
\end{align*}
$$
を用いる。したがって，Landau-Zenerモデルでは$ \varphi_s(\delta_{\mathrm{LZ}}) $，twisted Landau_Zenerモデルでは$ \varphi_s(\delta_{\mathrm{TLZ}})$と表せる。

## グラフ ($\kappa$を変化させる)
<img src="/resources/fig_es_Sp_k_0.01.png" alt="例">
<img src="/resources/fig_es_Sp_k_0.03.png" alt="例">
<img src="/resources/fig_es_Sp_k_0.05.png" alt="例">
<img src="/resources/fig_es_Sp_k_0.06.png" alt="例">
<img src="/resources/fig_es_Sp_k_0.07.png" alt="例">
<img src="/resources/fig_es_Sp_k_0.08.png" alt="例">
<img src="/resources/fig_es_Sp_k_0.09.png" alt="例">
<img src="/resources/fig_es_Sp_k_0.1.png" alt="例">
<img src="/resources/fig_es_Sp_k_0.2.png" alt="例">