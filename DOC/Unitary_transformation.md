# ユニタリ変換で$`x`$と$`z`$が入れ替え可能であることの証明
$`x`$成分と$`z`$成分はユニタリ変換によって入れ替えることが可能です。
このことは，遷移確率がHamiltonianの$`x \leftrightarrow z`$交換に依存しないことを意味します。

[!NOTE]
以下では，ユニタリ行列による相似変換のことをユニタリ変換と呼びます。    

## 結論
```math
U = 
\frac{1}{\sqrt{2}}
\begin{pmatrix}
    1 & 1\\
    1 & -1
\end{pmatrix}
```
によって$`x`$と$`z`$を入れ替えられます。

## 証明
### $`U`$がユニタリ行列であることの確認
```math
U U^{\dagger}
= \frac{1}{\sqrt{2}}
\begin{pmatrix}
    1 & 1\\
    1 & -1
\end{pmatrix}
\cdot
\begin{pmatrix}
    1 & 1\\
    1 & -1
\end{pmatrix}
= \frac{1}{2}
\begin{pmatrix}
    2 & 0\\
    0 & 2
\end{pmatrix}
=
\begin{pmatrix}
    1 & 0\\
    0 & 1
\end{pmatrix}
= I
```
また，$`U = U^{\dagger}`$より，$`U^{\dagger} U = I`$も明らかです。

### $`U`$によるユニタリ変換でPauli行列が入れ替わることの確認
#### $`x \rightarrow z`$
```math
\begin{align*}
    U^{\dagger} \sigma_x U
    &= \frac{1}{\sqrt{2}}
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}
    \cdot
    \begin{pmatrix}
        0 & 1\\
        1 & 0
    \end{pmatrix}
    \cdot
    \frac{1}{\sqrt{2}}
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}\\
    &= \frac{1}{2}
    \begin{pmatrix}
        1 & 1\\
        -1 & 1
    \end{pmatrix}
    \cdot
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}\\
    &=
    \frac{1}{2}
    \begin{pmatrix}
        2 & 0\\
        0 & -2
    \end{pmatrix}\\
    &=
    \begin{pmatrix}
        1 & 0\\
        0 & -1
    \end{pmatrix}\\
    &= \sigma_z
\end{align*}
```

#### $`y \rightarrow -y`$
```math
\begin{align*}
    U^{\dagger} \sigma_y U
    &= \frac{1}{\sqrt{2}}
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}
    \cdot
    \begin{pmatrix}
        0 & -i\\
        i & 0
    \end{pmatrix}
    \cdot
    \frac{1}{\sqrt{2}}
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}\\
    &= \frac{1}{2}
    \begin{pmatrix}
        i & -i\\
        -i & i
    \end{pmatrix}
    \cdot
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}\\
    &=
    \frac{1}{2}
    \begin{pmatrix}
        0 & 2i\\
        -2i & 0
    \end{pmatrix}\\
    &=
    \begin{pmatrix}
        0 & i\\
        -i & 0
    \end{pmatrix}\\
    &= -\sigma_y
\end{align*}
```

#### $`z \rightarrow x`$
```math
\begin{align*}
    U^{\dagger} \sigma_z U
    &= \frac{1}{\sqrt{2}}
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}
    \cdot
    \begin{pmatrix}
        1 & 0\\
        0 & -1
    \end{pmatrix}
    \cdot
    \frac{1}{\sqrt{2}}
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}\\
    &= \frac{1}{2}
    \begin{pmatrix}
        1 & -1\\
        1 & 1
    \end{pmatrix}
    \cdot
    \begin{pmatrix}
        1 & 1\\
        1 & -1
    \end{pmatrix}\\
    &=
    \frac{1}{2}
    \begin{pmatrix}
        0 & 2\\
        2 & 0
    \end{pmatrix}\\
    &=
    \begin{pmatrix}
        0 & 1\\
        1 & 0
    \end{pmatrix}\\
    &= \sigma_x
\end{align*}
```

### Hamiltonianの変換
```math
H = X \sigma_x + Y \sigma_y + Z \sigma_z
```
をユニタリ変換すると，確かに
```math
U^{\dagger} H U = X \sigma_z - Y \sigma_y + Z \sigma_x
```
となります。
Hamiltonianをユニタリ変換しても遷移確率は変わりません。
