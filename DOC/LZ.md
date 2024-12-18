# Landau-Zenerモデル
Hamiltonian
```math
H
=
\begin{pmatrix}
-\frac{1}{2} \nu t & -\frac{1}{2} \Delta \\
-\frac{1}{2} \Delta & \frac{1}{2} \nu t\\
\end{pmatrix}
```
において，upper stateからlower stateへ遷移する確率は，
```math
\mathcal{P}
= \exp(-2 \pi \delta)
```
で与えられます (Landau-Zener公式)。

## 他のモデルとの関係
### 本モデルを一般化したモデル
- [twisted Landau-Zenerモデル](/DOC/TLZ.md)
  - $`y`$成分を追加したモデル
- [multiple-passage Landau-Zenerモデル](/DOC/mpLZ.md)
  - LZモデルを複数回繰り返すモデル

## 証明方法
Landau-Zener公式は，1932年に，Majorana，Landau，Zener，Stückelbergによって独立に導出されました。

### 各方法の特徴
|証明者|方法|位相の計算|詳細|
|-|-|-|-|
|Majorana|Laplace変換|可|[こちら](/DOC/Majorana.md)|
|Landau|Dykhne公式|不可|coming soon|
|Zener|Weberの微分方程式，放物柱関数|可|[こちら](/DOC/Zener.md)|
|Stückelberg|WKB近似|不可|coming soon|
