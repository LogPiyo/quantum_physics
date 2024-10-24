# 目的
遷移確率を断熱パラメータの1次まで近似します。

# 計算

```math
\begin{align*}
P
&= \exp \left(-\frac{4}{|F|} \mathrm{Im} \, \int_0^{q_j} E(q) dq\right) \\
&= \exp \left(-\frac{4}{|F|} \mathrm{Im} \, \int_0^{q_j} \sqrt{X^2 + Y^2 + (Z - F \dot{\phi} / 2)^2} dq\right) \\
&\approx \exp \left(-\frac{4}{|F|} \mathrm{Im} \, \int_0^{q_j} \sqrt{E^2 -  Z F \dot{\phi}} dq\right) \\
&= \exp \left(-\frac{4}{|F|} \mathrm{Im} \, \int_0^{q_j} |E| \sqrt{1 -  \frac{Z F \dot{\phi}}{E^2}} dq\right) \\
&\approx \exp \left(-\frac{4}{|F|} \mathrm{Im} \, \int_0^{q_j} |E| (1 -  \frac{Z F \dot{\phi}}{2 E^2}) dq\right) \\
&= \exp \left(-\frac{4}{|F|} \mathrm{Im} \, \int_0^{q_j} |E| dq +  \frac{2}{|F|} \mathrm{Im} \, \int_0^{q_j} \frac{|E| Z F \dot{\phi}}{E^2} dq\right) \\
&= \exp \left(-\frac{4 (-F)}{|F|} \mathrm{Im} \, \int_0^{t_j} |E| dt +  \frac{2 (-F^2)}{|F|} \mathrm{Im} \, \int_0^{t_j} \frac{|E| Z \dot{\phi}}{E^2} dt \right) \\
\end{align*}
```