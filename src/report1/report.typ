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
ソースコード `src1.py` に基づいてプロットすると @figs:prob1-1 のようになる。

#figure(
  image("../../figs/prob1-1.png", width: 50%),
  caption: [$cos x$ の様々な有限差分法での微分比較],
) <figs:prob1-1>

これらの誤差のみを比較すると、@figs:prob1-2-1 のようになる。

#figure(
  image("../../figs/prob1-2-1.png", width: 50%),
  caption: [誤差の比較],
) <figs:prob1-2-1>

(2)
次に、ソースコード `src2.py` に基づいて厳密解と数値解の差をグリッド幅 $h$ の関数としてプロットすると @figs:prob1-2-2 のようになる。

#figure(
  image("../../figs/prob1-2-2.png", width: 50%),
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

== 課題3
