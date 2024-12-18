# Stokes位相のフィッティングで複数の候補が存在する問題について

<img src="/resources/double_passage_ex.png" alt="例" width=400>

上図で，解析解 (緑実線) と数値解 (青実線) が一致するようにStokes位相を決定したい。そこで，energy slope $\nu$ に対するStokes位相を示す，下図のようなグラフを出力する。
<br>
<img src="/resources/Fig_ex_es-Sp.png" alt="例" width=400>

このとき，1つのenergy slopeに対して複数の候補点が存在する。
そのため，現在は1つ前の点を基準として，連続につながるような1点を選択している。

## 複数の候補がある理由
### 理由1
```math
P
\rightarrow 4 \exp \left(-4 \frac{\omega}{|\omega|} \mathrm{Im} \int_0^{t_1} dt E'_+  \right) \left( 1 - \exp \left( -4 \frac{\omega}{|\omega|} \mathrm{Im} \int_0^{t_1} dt E'_+  \right) \right) \cos^2 \left(-\frac{\omega}{|\omega|} \int_{t_4}^{t_1} dt E'_+ + \varphi_s \right)
```
で，$\cos^2(\theta + n \pi) = \cos^2(\theta)$の周期性から，$\pi$ごとにフィッティングの候補点が存在する。この事実は，Stokes位相$\varphi_s$の各値に対して，

$$
\text{error} = \text{数値解} - \text{解析解}
$$
をプロットすると確かめられる。
<br>
例) $\nu = -80$の場合
<img src="/resources/fig_Sp-ae_v_-80.png" alt="例" width=400>
これは，Stokes位相がどんな値を取っても，緑実線 (解析解) が青実線 (数値解)を追い抜かないことを意味する。
### 理由2
理由1の例では極小値がちょうど0になっているが，パラメータを変えると極小値が0を下回ることがある。つまり，緑実線 (解析解) が青実線 (数値解)を追い越すようなパラメータ $\nu$ が存在する。
<br>
例1) $\nu = -50$の場合
<img src="/resources/fig_Sp-ae_v_-50.png" alt="例" width=400>
このパラメータのときは，$\pi$の移動に対して2点の候補が存在する。

例2) $\nu = -40$の場合
<img src="/resources/fig_Sp-ae_v_-40.png" alt="例" width=400>
緑実線 (解析解) と青実線 (数値解) の上下関係によってerrorの符号は異なる

## 対策
理由1: Stokes位相の探索幅は$\pi$より小さくする。
理由2: 2点の候補の平均を採用する