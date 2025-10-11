# 生物学の事例

到達目標
- 指数的増加と飽和効果をもつ成長モデルを理解する．
- 平衡点の解析を通じて安定性の概念を理解する．

## バクテリア増殖

### 実験的事実

- 栄養が十分にある培地に1匹のバクテリアを入れると，一定時間ごとに2分裂を繰り返して増殖する．
- 時間とともに個体数 $N(t)$ が急速に増加し，最初は**指数関数的に増加する**．
- しかし，資源の枯渇や老廃物の蓄積によって，最終的に増殖は**飽和**して一定値に近づく．

![動画]()

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

これを<span style="color:red">ロジスティックモデル</span>と呼ぶ．

## 微分方程式の性質

### 近似や極限を用いた解析

微分方程式をそのまま解析する前に，特別な場合や極端な場合での解の挙動が直観に合っているかを確認する．

| 条件 | 近似・極限 | 結論 |
| -- | -- | -- |
| $N \ll N_{\text{max}}$ | $(1 - \frac{N}{N_{\text{max}}}) \approx 1$ | 指数関数的増加 |
| $N \to N_{\text{max}}$ | $(1 - \frac{N}{N_{\text{max}}}) \to 0$ | 成長が止まる |

### 解

微分方程式

$$
\frac{dN}{dt} = rN\left(1 - \frac{N}{N_{\text{max}}}\right)
$$

変数分離して解く：

分母を部分分数分解：

$$
\frac{1}{\left(1 - \frac{N}{N_{\text{max}}}\right) N} = \frac{1}{N} + \frac{1/K}{1 - \frac{N}{N_{\text{max}}}}
$$

積分して：

$$
\int \left(\frac{1}{N} + \frac{1/K}{1 - \frac{N}{N_{\text{max}}}}\right)dN = \int rdt
$$

$$
\left(\log|N(t)| - \log|N_{\text{max}} - N(t)|\right) - \left(\log|N_0| - \log|N_{\text{max}} - N_0|\right) = rt.
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
t = np.linspace(0, 500, 200)
N = K / (1 + ((K - N0)/N0) * np.exp(-r * t))

fig, ax = plt.subplots(figsize=(5,3))
ax.plot(t, N, label="ロジスティック成長")
ax.axhline(K, color='gray', linestyle='--', label="平衡点 $N=K$")
ax.set_xlabel("時間 [h]")
ax.set_ylabel("個体数 $N(t)$")
ax.legend()
ax.grid(True)
fig.savefig("./3_logistic_growth.png", dpi=300)
plt.show()
```

```{note}
**演習1**

1. 次のパラメータを用いて，時刻 $t=0$ から $t=500$ までのロジスティック成長曲線を描け．
    - `r = 0.01`
    - `Nmax = 500`
    - `N_0 = 5`

2. `r` の値を2倍，`Nmax` の値を半分にしたとき，成長曲線の形はどう変化するか．
```

### 平衡点

平衡点（steady state）とは，成長が止まる点，すなわち

$$
\frac{dN}{dt} = 0
$$

を満たす $N$ を指す．
これは，式(1)より

$$
r(1 - \frac{N}{N_{\text{max}}})N = 0,
$$

であるから，これを満たす $N$ は

$$
N = 0 \quad \text{または} \quad N = N_{\text{max}}
$$

の2つ存在する．

### 平衡点の安定性解析

平衡点の安定性を調べるには，微分方程式の右辺の符号を調べる：

$$
f(N) = rN(1 - \frac{N}{N_{\text{max}}})
$$

| 区間          | (f(N)) の符号 | (N(t)) の変化 |
| ----------- | ---------- | ---------- |
| (N < 0)     | (f(N) > 0) | 増加（モデル外）   |
| (0 < N < K) | (f(N) > 0) | (N) は増加    |
| (N > K)     | (f(N) < 0) | (N) は減少    |

→ 図で示すと，(N=K) に向かって矢印が集まる．

したがって：

- (N=0)：**不安定平衡点**（少し増えると増殖へ）
- (N=K)：**安定平衡点**（収束点，飽和状態）

```{note}
**演習2**

1. 平衡点 $N=0$ と $N=N_{\text{max}}$ の安定性を，グラフ上の矢印や $f(N)$ の符号を使って説明せよ．
```

---

## 7. まとめ

| 項目    | 内容                                                                           |
| ----- | ---------------------------------------------------------------------------- |
| モデル   | (\displaystyle \frac{dN}{dt} = rN\left(1 - \frac{N}{K}\right))               |
| パラメータ | (r)：成長率，(K)：環境収容力                                                            |
| 平衡点   | (N=0)（不安定），(N=K)（安定）                                                         |
| 解     | (\displaystyle N(t) = \frac{K}{1 + \left(\frac{K - N_0}{N_0}\right)e^{-rt}}) |
| 振る舞い  | 初期は指数増加，後に飽和（S字型曲線）                                                          |
| 意義    | 生物集団の増殖・感染拡大・資源利用などの基本モデル                                                    |

---

### 💡補足：S字カーブの特徴

- 初期：資源が豊富 → 指数的増加
- 中期：競争・制限効果 → 増加が緩やか
- 後期：飽和 → 定常状態 (N=K)

---


## まとめ

本日作成したipynbファイルに追記する形で次の課題に回答し，WebClassの第3回課題から提出せよ．

```{note}
**復習**

(1) 本日学んだことを3つ箇条書きで述べよ．
```

### 参考文献

- P.F. Verhulst, *Notice sur la loi que la population poursuit dans son accroissement*, 1838.
- M. Murray, *Mathematical Biology*, Springer.
- 伊藤清三『数理生物学入門』東京大学出版会.
