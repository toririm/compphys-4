#set text(font: "Hiragino Mincho ProN")
#let mixed(body) = {
  set text(weight: "extrabold")
  set text(font: "Hiragino Sans", weight: 600)
  body
}

#show heading: mixed

= 計算物理学Ⅳ レポート2
202210807 時長隆乃介

== 速度Verlet法の導出
$v(t + Delta t)$ を $t$ の周りで Taylor 展開すると
$
  v(t + Delta t)
  &= v(t)
  + (d v(t))/(d t) Delta t
  + 1/2 (d^2 v(t))/(d t^2) Delta t^2
  + O(Delta t^3) \

  &= v(t)
  + 1/m F(t) Delta t
  + 1/(2m) (d F(t))/(d t) Delta t^2
  + O(Delta t^3) \

  &= v(t)
  + 1/m F(t) Delta t
  + 1/(2m) (F(t + Delta t) - F(t)) Delta t
  + O(Delta t^3) \

  &= v(t)
  + 1/(2m) (F(t) + F(t + Delta t)) Delta t
  + O(Delta t^3) \
$
よって導出された。

== Morseポテンシャルから力の導出
Morseポテンシャル
$
  V(r) = D_e (1 - e^(-a(r - r_e)))^2
$
を $r$ で微分して
$
  F(r) &= -(d V(r))/(d r) \
  &= 2 D_e (1 - e^(-a(r - r_e))) dot (a e^(-a(r-r_e)))
$

== `LJ.c`の実験
`LJ.c` を 初期温度 0-1000K の間で 100K 区切りで実行し、その結果を元に拡散係数を求め、プロットした。
グラフを @fig:plot に示す。

#figure(
  image("plot.png", width: 80%),
  caption: [拡散係数の温度依存性]
) <fig:plot>
