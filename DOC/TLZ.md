# twisted Landau-Zenerモデル

Hamiltonian
```math
\begin{align*}
     H_\mathrm{TLZ}(t)
     &=
     \begin{pmatrix}
          m & \nu F t - i \frac{1}{2} \kappa_{||} \nu^2 (F t)^2\\
          \nu F t + i \frac{1}{2} \kappa_{||} \nu^2 (F t)^2 & m
     \end{pmatrix} \\
     &=
     \begin{pmatrix}
          \Delta_z & \varepsilon_0 \omega t - i \frac{1}{2} \Delta_y (\omega t)^2\\
          \varepsilon_0 \omega t + i \frac{1}{2} \Delta_y (\omega t)^2 & \Delta_z
     \end{pmatrix} 
\end{align*}
```
を**twisted Landau-Zenerモデル**と呼びます。このモデルのupper stateからlower stateへ遷移する確率は，
```math
P(F)
= \exp \left[-\pi \frac{(\Delta_0 - \frac{\kappa_g \varepsilon_0 F}{4})^2}{\varepsilon_0 |F|}\right]
```
で与えられます。ここで，
```math
\Delta_0 - \frac{\kappa_g \varepsilon_0 F}{4} = 0,
```
すなわち，遷移確率が1となるような有限の$`F`$が存在することが，Landau-Zenerモデルとの違いです。この現象を**完全トンネル** (perfect tunneling) と呼びます。

## 遷移確率の証明法