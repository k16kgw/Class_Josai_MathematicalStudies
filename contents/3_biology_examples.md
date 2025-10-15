# 生物学の事例

到達目標
- 指数的増加と飽和効果をもつ成長モデルを理解する．
- 平衡点の解析を通じて安定性の概念を理解する．

キーワード
- 飽和効果
- 安定性

準備
1. anacondaを使用し，<span style="color:red">jupyter lab</span>を起動する．
2. `Documents（書類）/mathematical_studies`フォルダをダブルクリックで開き`+`をクリックして新しいファイルを作成する．
3. ファイル名を`3_{学籍番号}_{氏名}.ipynb`に変更する．例：`3_SI25999_香川渓一郎.ipynb`


## バクテリア増殖

### 実験的事実

- 栄養が十分にある培地に1匹のバクテリアを入れると，一定時間ごとに分裂を繰り返して増殖する．
- 時間とともに個体数 $N(t)$ が急速に増加し，最初は**指数関数的に増加する**．
- しかし，資源の枯渇や老廃物の蓄積によって，最終的に増殖は**飽和**して一定値に近づく．

<!-- ![バクテリアコロニーの成長を撮影した動画 [1]](/contents/figs/3/elife-48885-video1.mp4) -->

<iframe width="560" height="315" src="https://www.youtube.com/embed/MUl1R_tho5g?si=mVkTBQoLgVBPKtSS" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

[1]より

<iframe width="560" height="315" src="https://www.youtube.com/embed/z43HNp3zVAA?si=KnnaOC3PeW-HIcpz" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

[2]より

### 指数的増加モデルの復習

資源が無限にあると仮定した場合，時刻 $t$ の個体数を $N(t)$，増殖率を $r > 0$ として個体数の時間変化は次の式で表される．

$$
\frac{dN}{dt} = rN.
$$

初期時刻 $t=0$ での個体数を $N_0$ とすれば，この微分方程式の解は次のように指数関数で表される．

$$
N(t) = N_0 e^{rt}
$$

つまり，個体数は初期値から無限に増加していく．
しかし<span style="color:red">現実には増殖の限界がある</span>．

### 飽和効果を考慮したモデル化

- 個体数が大きくなると，栄養や空間が不足して増殖が抑制される．
- 増殖率が個体数の増加とともに減少すると考える．
- 例えば個体数が環境に収容できる上限の値 $N_{\text{max}} > 0$ に届いたら増殖率を $0$ にする．

$$
\frac{dN}{dt} = r \Bigl(1 - \frac{N}{N_{\text{max}}}\Bigr) N.
\tag{1}
$$

この数理モデルを<span style="color:red">ロジスティックモデル</span>と呼び，微分方程式 (1) を<span style="color:red">ロジスティック方程式</span>と呼ぶ．

## 微分方程式の性質

### 近似や極限を用いた解析

微分方程式をそのまま解析する前に，特別な場合や極端な場合での解の挙動が直観に合っているかを確認する．

| 条件 | 近似・極限 | 微分方程式 | 結論 |
| -- | -- | -- | -- |
| $N \ll N_{\text{max}}$ | $(1 - \frac{N}{N_{\text{max}}}) \approx 1$ | $\frac{dN}{dt} = rN$ | 指数関数的増加 |
| $N \to N_{\text{max}}$ | $(1 - \frac{N}{N_{\text{max}}}) \to 0$     | $\frac{dN}{dt} = 0$  | 成長が止まる |

### 微分方程式の解

ロジスティック方程式

$$
\frac{dN}{dt} = rN\left(1 - \frac{N}{N_{\text{max}}}\right)
$$

を変数分離して解く．
両辺を$N\left(1 - \frac{N}{N_{\text{max}}}\right)$で除し，
$\frac{1}{N\left(1 - \frac{N}{N_{\text{max}}}\right)}$を部分分数分解すると

$$
\frac{1}{N\left(1 - \frac{N}{N_{\text{max}}}\right)} = \frac{1}{N} + \frac{\frac{1}{N_{\text{max}}}}{1 - \frac{N}{N_{\text{max}}}}
$$

であるから，微分方程式は

$$
\left( \frac{1}{N} + \frac{\frac{1}{N_{\text{max}}}}{1 - \frac{N}{N_{\text{max}}}} \right) \frac{dN}{dt} = r.
$$

これを時間に関して区間$(0,t)$で積分すれば

$$
\int_0^t \left(\frac{1}{N} + \frac{\frac{1}{N_{\text{max}}}}{1 - \frac{N}{N_{\text{max}}}}\right)dN = \int_0^t rdt
$$

$$
\Bigl(\log|N(t)| - \log|1 - \frac{N(t)}{N_{\text{max}}}|\Bigr) - \Bigl(\log|N_0| - \log|1 - \frac{N_0}{N_{\text{max}}}|\Bigr) = rt.
$$

従って

$$
N(t) = \frac{N_{\text{max}}}{1 + \left(\frac{N_{\text{max}} - N_0}{N_0}\right)e^{-rt}}.
\tag{2}
$$

### グラフのプロット

準備
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'
```

変数設定
```python
r = 0.02    # 成長率
Nmax = 1000    # 環境収容力
N0 = 10     # 初期個体数
```

プロット
```python
t = np.linspace(0, 800, 200)
N = Nmax / (1 + ((Nmax - N0)/N0) * np.exp(-r * t))

fig, ax = plt.subplots(figsize=(5,3))
ax.plot(t, N)
ax.axhline(Nmax, color='gray', label=f"$N=Nmax={Nmax}$") # Nmaxのラインを描画
ax.set_xlabel(r"時刻 $t$")
ax.set_ylabel(r"個体数 $N(t)$")
ax.legend() # labelの表示
ax.grid(True) # 格子の表示
plt.tight_layout() # レイアウトを調整
fig.savefig("./3_logistic_growth.png")
plt.show()
```

![グラフのプロット](/contents/figs/3/logistic_growth.png)

ここで示した曲線を<span style="color:red">ロジスティック曲線</span>と呼ぶ．

```{note}
**演習1**

1. 次のパラメタの下でロジスティック成長曲線を描け．ただし時刻は$t=0$から成長が停止したと思われる時刻までを自分で設定せよ．
    - `r = 0.01`
    - `Nmax = 500`
    - `N0 = 5`

2. `r` の値を2倍，`Nmax` の値を半分にしたとき，成長曲線の形はどう変化するか．グラフを2つ重ねて描画して比較せよ．
```

### 平衡点

平衡点（steady state）：成長が止まる点，すなわち$\frac{dN}{dt} = 0$を満たす $N=N_{\text{eq}}$．

ロジスティック方程式 (1) の平衡点は

$$
r\left(1 - \frac{N}{N_{\text{max}}}\right)N = 0
$$

を満たす $N$ であるから，

$$
N_{\text{eq}} = 0 \quad \text{または} \quad N = N_{\text{max}}
$$

の2つである．

### 平衡点の安定性解析

微分方程式の平衡点は，**安定性**の観点から2種類に分類できる．

- <span style="color:red">安定平衡点</span>：そこから僅かにズレた値を初期値とすると，解は安定平衡点に限りなく近づく．
- <span style="color:red">不安定平衡点</span>：そこから僅かにズレた値を初期値とすると，解は不安定平衡点から離れていく．

平衡点の安定性を調べるには，微分方程式の右辺の符号を調べれば良い．
右辺を次のように$f(N)$と置く．

$$
f(N) = rN(1 - \frac{N}{N_{\text{max}}})
$$

| 区間                      | $f(N)$の符号 | (N(t)) の変化 |
| -----------              | ----------  | ----------   |
| $N < 0$                  | $f(N) > 0$  | 増加（モデル外）|
| $0 < N < N_{\text{max}}$ | $f(N) > 0$  | $N$は増加     |
| $N_{\text{max}} < N$     | $f(N) < 0$  | $N$は減少     |

→ 増減の方向を図で示すと，$N=0$からは離れる方向に$N$は変化し，$N=N_{\text{max}}$に近づく方向に$N$は変化することが分かる．

従って

- $N=0$：**不安定平衡点**（少し増えると増殖へ）
- $N=N_{\text{max}}$：**安定平衡点**（少し減らしても収束する．飽和状態）

### グラフのプロット

複数の値の初期値$N_0$についてロジスティック曲線をプロットし，曲線が $N=0$ から離れ，$N=N_{\text{max}}$ に漸近することを確認する．

変数の設定
```python
r = 0.02    # 成長率
Nmax = 1000    # 環境収容力
N0_list = [0, 1, 10, 100, 500, 900, 990, 999, 1000]
```

グラフのプロット
```python

t = np.linspace(0, 800, 200)

fig, ax = plt.subplots(figsize=(5,3))
ax.axhline(Nmax, color='gray') # Nmaxのラインを描画

for N0 in N0_list:
    N = Nmax*N0 / (N0 + (Nmax - N0) * np.exp(-r * t))
    # ax.plot(t, N, label=f"$N0={N0}$", color=cmap(norm(N0)))
    ax.plot(t, N, label=f"$N0={N0}$")
ax.set_xlabel(r"時刻 $t$")
ax.set_ylabel(r"個体数 $N(t)$")
ax.legend() # labelの表示
ax.grid(True) # 格子の表示
plt.tight_layout() # レイアウトを調整
fig.savefig("./3_multi_logistic_growth.png")
plt.show()
```

![複数の初期値に対するロジスティック曲線](/contents/figs/3/multi_logistic_growth.png)

```{note}
**演習2**

次のパラメタの下で，初期値$N_0$を複数の値でプロットし，平衡点 $N=0$ と $N=N_{\text{max}}$ の安定性を曲線の形から確認せよ，
    - `r = 0.01`
    - `Nmax = 500`
```

## タンチョウの個体数の変化をロジスティックモデルで説明する

![タンチョウの個体数の変化](/contents/figs/3/tancho.png)

## まとめ

| 項目      | 内容                                                                          |
| -----    | ---------------------------------------------------------------------------- |
| モデル    | $\displaystyle \frac{dN}{dt} = rN\left(1 - \frac{N}{N_{\text{max}}}\right)$  |
| パラメタ  | $r$：成長率，$N_{\text{max}}$：環境収容力                                        |
| 平衡点    | $N=0$（不安定），$N=N_{\text{max}}$（安定）                                     |
| 解       | $\displaystyle N(t) = \frac{N_{\text{max}}}{1 + \left(\frac{N_{\text{max}} - N_0}{N_0}\right)e^{-rt}}$ |
| 振る舞い  | 初期は指数関数的増加，後に飽和し，S字型曲線をなす．                                   |
| 意義     | 生物集団の増殖・感染拡大・資源利用などの基本モデル                                    |

### 提出課題

1. 本日作成したipynbファイルをWebClassの「第3回課題」から提出せよ．

## 参考文献

[1] Liyang Xiong, Yuansheng Cao, Robert Cooper, Wouter-Jan Rappel, Jeff Hasty, Lev Tsimring, "Flower-like patterns in multi-species bacterial colonies," eLife, 9 (2020) e48885.

[2] HHMI BioInteractive Video, "Bacterial Growth," URL:https://www.biointeractive.org/classroom-resources/bacterial-growth （閲覧日：2025.10.15）

[3] 小寺隆幸, 数学で考える環境間題, 明治図書 (2004).

[4] 梅野善雄, "タンチョウの個体数変化とロジスティック曲線," 数学教育学会誌, 50 (2009) 5–13. https://doi.org/10.34323/mesj.50.1-2_5.
  

