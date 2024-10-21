#set text(font: "Noto Serif JP")
#let mixed(body) = {
  set text(weight: "extrabold")
  set text(font: "Noto Sans JP", weight: "bold")
  body
}

#show heading: mixed

= 計算物理学IV レポート課題1

時長隆乃介 202210807

== 課題1

(1) 
ソースコード #link("https://google.com")[`src1.py`] に基づいてプロットすると @figs:prob1-1 のようになる。

#figure(
  image("figs/prob1-1.png", width: 50%),
  caption: [$cos x$ の様々な有限差分法での微分比較],
) <figs:prob1-1>

これらの誤差のみを比較すると、@figs:prob1-2-1 のようになる。

#figure(
  image("figs/prob1-2-1.png", width: 50%),
  caption: [誤差の比較],
) <figs:prob1-2-1>

#pagebreak()

(2)
次に、ソースコード `src2.py` に基づいて厳密解と数値解の差をグリッド幅 $h$ の関数としてプロットすると @figs:prob1-2-2 のようになる。

#figure(
  image("figs/prob1-2-2.png", width: 50%),
  caption: [厳密解と数値解の差],
) <figs:prob1-2-2>

このプログラムを用いて、平均残差が $10^(-6)$ 以下になるときの $Delta x$ を見積もると @tab:prob1-2-3 のようになる。

#figure(
  table(
    columns: 3,
    table.header([手法], [$Delta x$], [$log_10 (Delta x)$]),
    table.hline(),
    [前方差分], [$2.422 times 10^(-6)$], [$-5.616$],
    [後方差分], [$2.422 times 10^(-6)$], [$-5.616$],
    [中央差分], [$3.075 times 10^(-3)$], [$-2.512$]
  ),
  caption: [平均残差が $10^(-6)$ となる $Delta x$ の見積もり],
) <tab:prob1-2-3>


== 課題2
$
  (partial^2 f(x))/(partial x^2) = (partial f'(x))/ (partial x)
$
$f'(x)$ に対する中央差分を用いて
$
  (partial f'(x))/ (partial x) = (f'(x + h) - f'(x - h)) / (2h)
$
$f'(x + h)$ に対して前方差分、$f'(x - h)$ に対して後方差分を用いると
$
  (f'(x + h) - f'(x - h)) / (2h) = ({f(x + 2h) - f(x)}slash 2h - {f(x) - f(x - 2h)}slash 2h) / (2h)
$
$2h$ を $h$ として置きなおすと
$
  (partial^2 f(x))/(partial x^2) = (f(x + h) - 2f(x) + f(x - h)) / h^2
$

#pagebreak()

== 課題3
(1)
有限区間 $[0, 2pi]$ で考える。 
$V(x) = 0$ として与式を書き下すと、一次元の時間に依存しないシュレディンガー方程式は
$
  -1/2(partial^2 psi(x))/(partial x^2) = E psi(x)
$
=== 厳密解
一般解は
$
  psi(x) = A cos(k x) + B sin(k x)
$
で与えられる。
与式に代入して
$
  1/2 k^2 = E space therefore k = plus.minus sqrt(2E)
$
境界条件 $psi(0)=psi(2pi)=0$ より
$
  psi(0) = A = 0 \
  psi(2pi) = A cos(2pi k) + B sin(2pi k) = 0
$
したがって、非自明な解として
$
  sin (2pi k) = 0, space 2pi k = n pi \
  k = n/2, space therefore E = n^2 / 8 space (n = 1, 2, 3, ...)
$
対応する波動関数は
$
  psi_n (x) = B sin((n x)/2)
$

規格化条件より
$
  B = 1 / sqrt(pi), space therefore psi_n (x) = 1 / sqrt(pi) sin((n x)/2)
$

=== 数値解
有限差分法を用いると、与式は
$
  -1/2(psi(x + h) - 2 psi(x) + psi(x - h)) / h^2 = E psi(x)
$
ハミルトニアン行列 $H$ は
$
  H_(i, j) = -1/2 (delta_(i, j + 1) - 2 delta_(i, j) - 1/2 delta_(i, j - 1))/h^2
$
として
$
  H mat(psi(x_1); psi(x_2); dots.v; psi(x_N)) = E mat(psi(x_1); psi(x_2); dots.v; psi(x_N))
$
よって、行列 $H$ の 固有値 $E$ と固有ベクトル $psi$ を求めることで数値解を得ることができる。

`src3.py`によって求めて、厳密解と共にプロットしたものが @figs:prob1-3 である。

#figure(
  image("figs/prob1-3.png", width: 50%),
  caption: [厳密解と数値解の比較],
) <figs:prob1-3>

なぜか第二励起状態では符号が反転している。たぶんeigenvectorsの取り方が悪いのだと思う。

固有値は次の出力を得た。
```
numerical: [0.12450139 0.49800435 1.1205052 ]
analytical: [0.125, 0.5, 1.125]
```
