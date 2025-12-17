# BZ反応のモデル

到達目標
- BZ反応に代表されるリミットサイクル振動を理解する．
- 数値シミュレーションを活用して周期解を理解する．

準備
1. anacondaを使用し，<span style="color:red">jupyter lab</span>を起動する．
2. `Documents（書類）/mathematical_studies`フォルダをダブルクリックで開き`+`をクリックして新しいファイルを作成する．
3. ファイル名を`8_{学籍番号}_{氏名}.ipynb`に変更する．例：`8_SI25999_香川渓一郎.ipynb`

## BZ反応（Belousov–Zhabotinsky reaction）

- 同じ溶液の中にある化学反応が**周期的に色を変え続ける現象**．
- 1950年代に発見される．
- 多くの化学反応は「反応物が減り生成物が増えたら終わる」が，BZ反応は**外部入力なしでも勝手に振動**する．

### 実験動画

<iframe width="560" height="315" src="https://www.youtube.com/embed/Nrwc8Ojvnxk?si=cKTVk1W7Grgsnzs1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

[「振動反応　～BZ(ベローゾフ・ジャボチンスキー)反応～」VCPteam’s blog](https://vcpteam.hatenablog.com/entry/2022/07/05/215333)

### 代表的なモデル

- Field–Körös–Noyes(FKN)モデル
- オレゴネータ(Oregonator)モデル
- Keener–Tysonモデル

## FKNモデル

次の5つの化学成分が重要な役割を果たしている．

$$
& X = \mathrm{HBrO_2}, 
\quad Y = \mathrm{Br^{-}},
\quad Z = \mathrm{Ce^{4+}},
\\
& A = \mathrm{BrO_3^{-}},
\quad P = \mathrm{HOBr},
$$

これらの間に次の化学反応が生じている．

$$
& A + Y \xrightarrow{k_1} X + P,
\\
& X + Y \xrightarrow{k_2} 2P,
\\
& A + X \xrightarrow{k_3} 2X + 2Z,
\\
& 2X \xrightarrow{k_4} A + P,
\\
& Z \xrightarrow{k_5} f Y,
$$

ただし$k_i$($i=1,\ldots,5$)は反応速度定数，$f$は化学量論因子で$f=0.5$とされることが多い．

## Oregonatorモデル

FKNモデルを3つの成分に帰着したモデル．
Fieldらがオレゴン大学で研究していたことから名付けられた．

$X=\mathrm{HBrO_2}$，$Y=\mathrm{Br^{-}}$，$Z=\mathrm{Ce^{4+}}$の3成分に注目し，時刻$t$におけるそれぞれの濃度を次のようにおく．
- $x(t)$：中間体$\mathrm{HBrO_2}$の濃度
- $y(t)$：還元剤$\mathrm{Br^{-}}$の濃度
- $z(t)$：触媒$\mathrm{Ce^{4+}}$の濃度

また$A = \mathrm{BrO_3^{-}}$はほとんど定数$a$であると仮定し，$X,Y,Z$の3成分の時間変化を考えると

$$
\begin{cases}
\dfrac{dx}{dt} = k_1 ay - k_2 xy + k_3 ax - k_4 x^2,
\\
\dfrac{dy}{dt} = -k_1 ay - k_2 xy + f k_5 z,
\\
\dfrac{dz}{dt} = 2k_3 ax - k_5 z,
\end{cases}
$$

となる．
<!-- 
ここでは$z$がほとんど変化しない状況を考え，次の簡略版Oregonatorモデルを扱う．

$$
\begin{cases}
\dfrac{dx}{dt} = q y - x y + x(1 - x),
\\
\dfrac{dy}{dt} = - q y - x y + f,
\end{cases}
$$

ここに
- $q>0$：反応比
- $f>0$：供給（定数）
- 非線形項は$xy$, $x^2$の2つ
- 適切なパラメタのもとでは，静止せずに閉じた軌道（リミットサイクル）に沿って動く．

## 平衡点の解析

平衡点は$\dfrac{dx}{dt}=\dfrac{dy}{dt}=0$で求まるが，ここでは表式が複雑になるため数値的に求める．

```python
import mpmath as mp
q, f = 0.002, 1.1

# f1(x, y) = 0, f2(x, y) = 0 を解きたい
def f1(x, y):
    return q*y - x*y + x*(1 - x)
def f2(x, y):
    return -q*y - x*y + f
# 初期値 (x0, y0) = (0.2, 1.0) からニュートン法で解く
x_star, y_star = mp.findroot((f1, f2), (0.2, 1.0))

print("平衡点: x* =", x_star, ", y* =", y_star)
```
 -->
<!-- 
```python
import mpmath as mp
q, f = 0.002, 1.1

x_star, y_star = mp.findroot(
    (lambda x, y: q*y - x*y + x*(1-x),
     lambda x, y: -q*y - x*y + f),
    (0.2, 1.0)
)
print("平衡点: x* =", x_star, ", y* =", y_star)
```
 -->
<!-- 
この平衡点は近くを初期値とすると
- 内側からは外へ
- 外側からは内へ
と向かい，**ある閉じた軌道に吸い込まれていく**．
この軌道を<span style="color:red">リミットサイクル</span>と呼ぶ．

## 数値シミュレーション

右辺の実装
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'

q = 0.002
f = 1.1

def oregonator_rhs(x, y):
    dx = q*y - x*y + x*(1-x)
    dy = -q*y - x*y + f
    return dx, dy
```


```python
def euler_oregonator(f, x0, y0, t):
    X = np.empty_like(t); Y = np.empty_like(t)
    X[0], Y[0] = x0, y0
    h = t[1] - t[0]

    for k in range(len(t) - 1):
        dx, dy = f(X[k], Y[k])   # 現在値で傾きを評価
        X[k+1] = X[k] + h * dx   # x を更新
        Y[k+1] = Y[k] + h * dy   # y を更新

    return X, Y
```

計算の実施
```python
t = np.linspace(0, 200, 10000)  # h が小さめになるように刻みをとる
x0, y0 = 0.2, 1.2
X, Y = euler_oregonator(oregonator_rhs, x0, y0, t)
```

可視化（相平面）
```python
plt.figure(figsize=(5,4))
plt.plot(X, Y, 'b-')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Oregonator (BZ reaction) のリミットサイクル")
plt.grid(True, alpha=0.3)
plt.show()
```

- 内側の軌道は外側に発散
- 外側の軌道は内側へ収束
- 最終的に一定の閉じたループへ吸い込まれる：リミットサイクル（limit cycle）

時系列表示
```python
plt.figure(figsize=(6,3))
plt.plot(t, X, label="x(t)")
plt.plot(t, Y, label="y(t)")
plt.xlabel("t"); plt.ylabel("濃度")
plt.legend(); plt.grid(True, alpha=0.3)
plt.title("BZ反応の濃度振動")
plt.show()
```

結果：
- $x(t)$と$y(t)$は周期的に振動する．
- 化学反応が時間に関して周期的な挙動を示す現象に相当する．

## リミットサイクルの数学的特徴

- リミットサイクルは

  $$
  \gamma(t):\mathbb{R}\to\mathbb{R}^2,\quad \gamma(t+T)=\gamma(t)
  $$

  という周期解で，その軌道の近くにある点が時間とともに**吸い寄せられる**（安定極限周期）．
- BZ反応は非線形項とフィードバックによって，この構造を内包している．

構造
- 内側からくる軌道：外へ押される（不安定）
- 外側からくる軌道：内側へ引っ張られる（安定）
- 両者の境界として閉曲線（安定なリミットサイクル）が出現


```{note}
**演習1**

Oregonatorモデルを用いて異なる初期値から出発しても同じリミットサイクルに収束することを示せ．
```

```{note}
**演習2**

パラメタ$q$と$f$の値によっては安定なリミットサイクルが消滅することがある．パラメタ$q$と$f$を変化させ，安定なリミットサイクルが消滅・出現する様子を調べよ．
```
 -->

## Keener–Tysonモデル

特異摂動と呼ばれる手続きによって変数$y$の時間変化が分からないほど遅いとすることで，Oregonatorモデルは次のように簡略化（縮約）することができる．

$$
\begin{cases}
\varepsilon \dfrac{dx}{dt} = fz\dfrac{q-x}{q+x} + x(1-x),
\\
\dfrac{dz}{dt} = x - z,
\end{cases}
$$

この2変数系について数値シミュレーションによって解のダイナミクスを調べてみる．

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'

eps = 0.01
q = 0.01
f = 0.8

def kt_rhs(x, z):
    dx = (f*z*(q-x)/(q+x) + x*(1-x))/eps
    dz = x - z
    return dx, dz

def euler_kt(f, x0, z0, t):
    X = np.empty_like(t); Z = np.empty_like(t)
    X[0], Z[0] = x0, z0
    h = t[1] - t[0]

    for k in range(len(t) - 1):
        dx, dz = f(X[k], Z[k])   # 現在値で傾きを評価
        X[k+1] = X[k] + h * dx   # x を更新
        Z[k+1] = Z[k] + h * dz   # z を更新

    return X, Z
```

相図の作成
```python
x = np.linspace(0, 1.2, 25)
z = np.linspace(0, 0.4, 25)
X, Z = np.meshgrid(x, z)
U, V = kt_rhs(X, Z)

fig, ax = plt.subplots(figsize=(5,4))
ax.quiver(X, Z, U, V, color='gray', angles='xy')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("KTモデルの相図")
ax.grid(True, alpha=0.3)
plt.savefig(f"8_kt_phase_eps={eps}_q={q}_f={f}.png")
plt.show()
```

計算の実施と軌道のプロット
```python
t = np.linspace(0, 20, 1000000) # 分割幅が小さくなるよう注意
x0, z0 = 0.2, 1.2
X, Z = euler_kt(kt_rhs, x0, z0, t)

plt.figure(figsize=(5,4))
plt.plot(X, Z, 'b-')
plt.scatter(x0, z0, c='r')
plt.xlabel("x")
plt.ylabel("z")
plt.title("Keener–Tysonモデルのリミットサイクル")
plt.grid(True, alpha=0.3)
plt.savefig(f"8_kt_xz_orbit_eps={eps}_q={q}_f={f}.png")
plt.show()
```

上記コードを実行すると次のような図が得られる．

![リミットサイクル](/contents/report/8/8_kt_xz_orbit_eps=0.01_q=0.01_f=0.8_x0=0.2_z0=0.2.png)

- $xz$平面内で初期値から始まった解の軌道はある一定の閉じたループへ吸い込まれる様子が見られる．
- このループを<span style="color:red">リミットサイクル</span>(limit cycle)と呼ぶ．

```{note}
<span style="color:red">**課題1**</span>

上記，軌道をプロットするコードで，`x0`，`z0`の値を変えても同じループに吸い込まれることを確認せよ．

特にループの外側の点を初期値とした場合にループの外側から吸い込まれる図を出力せよ．
```

$x$と$z$の時間変化
```python
plt.figure(figsize=(6,3))
plt.plot(t, X, label="x(t)")
plt.plot(t, Z, label="z(t)")
plt.xlabel("t"); plt.ylabel("濃度")
plt.legend(); plt.grid(True, alpha=0.3)
plt.title("Keener–Tysonモデルの濃度振動")
plt.savefig(f"8_kt_t-xz_eps={eps}_q={q}_f={f}.png")
plt.show()
```

上記コードを実行すると次のような図が得られる．

![リミットサイクル](/contents/report/8/8_kt_t-xz_eps=0.01_q=0.01_f=0.8_x0=0.2_z0=0.2.png)

解の時間変化をプロットすると上の図のように周期的に振動していることが確認できる．
この周期的な挙動が，解がリミットサイクル上を回転していることに対応している．

### 課題

```{note}
<span style="color:red">**課題2**</span>

Keener–Tysonモデルでは`eps`，`q`，`r`の3つのパラメタが存在するが，このパラメタの値を変化させることで解の軌道がどのように変化するかを確かめるために，次の3つの問いに答えよ．

1. `eps`の値を$0.01$から$0.1, 1, 10$と変化させたときの軌道と時間変化がどのように変化するか，図と共に述べよ．ただし，`q=0.01, r=0.8`とする．
2. `q`の値を$0.03, 0.05, 0.061, 0.07$と変化させたときの軌道と時間変化がどのように変化するか，図と共に述べよ．ただし，`eps=0.01, r=0.8`とする．
3. `r`の値を自由に変化させ，リミットサイクルが得られる値の範囲を調べよ．`r`の値を動かすのは小数点以下第一位までで良い．ただし，`eps=0.01, q=0.01`とする．
```

---

## 応用

### 自律振動

外部入力なしでも周期運動が起きる
→ 生体リズム（心臓鼓動，概日リズム）の研究モデル

### 化学パターン形成

BZ反応に空間的な構造も含めて考えると
- 同心円波
- スパイラル波
- 迷路状パターン
が自発的に形成される．

[「研究内容」末松　J.　信彦](https://www.isc.meiji.ac.jp/~suematsu/research/pattern.html)

### 脳科学・計算科学への応用

- 神経活動のモデル化
- パターン形成と自己組織化の概念学習
- ケミカルコンピューティング（cf.[「物質を使った新しい計算原理に関する情報科学 (Physical and Chemical Computing & Molecular Intelligence)」瀧ノ上研究室](https://takinoue-lab.jp/research/research-computing/)）

---

## まとめ

| 項目         | 内容                      |
| ---------- | ----------------------- |
| BZ反応       | 化学的非平衡系で生じる自律振動         |
| Oregonator | BZ反応を説明する簡略モデル          |
| リミットサイクル   | 吸引性を持つ周期軌道              |

パラメタ$q$と$f$の値によっては安定なリミットサイクルが消滅することがある．パラメタ$q$と$f$を変化させ，安定なリミットサイクルが消滅・出現する様子を調べよ．
3. （発展）数値的に周期 (T) を推定せよ（ゼロクロス法など）。


