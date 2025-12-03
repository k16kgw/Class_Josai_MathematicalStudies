# BZ反応のモデル

到達目標
- BZ反応に代表されるリミットサイクル振動を理解する．
- 数値シミュレーションを活用して周期解を理解する．

準備
1. anacondaを使用し，<span style="color:red">jupyter lab</span>を起動する．
2. `Documents（書類）/mathematical_studies`フォルダをダブルクリックで開き`+`をクリックして新しいファイルを作成する．
3. ファイル名を`8_{学籍番号}_{氏名}.ipynb`に変更する．例：`3_SI25999_香川渓一郎.ipynb`

## BZ反応（Belousov–Zhabotinsky reaction）

- 同じ溶液の中にある化学反応が**周期的に色を変え続ける現象**．
- 1950年代に発見される．
- 多くの化学反応は「反応物が減り生成物が増えたら終わる」が，BZ反応は**外部入力なしでも勝手に振動**する．

### 実験動画

<iframe width="560" height="315" src="https://www.youtube.com/embed/Nrwc8Ojvnxk?si=cKTVk1W7Grgsnzs1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

[「振動反応　～BZ(ベローゾフ・ジャボチンスキー)反応～」VCPteam’s blog](https://vcpteam.hatenablog.com/entry/2022/07/05/215333)

## 化学反応系と非線形性

化学反応系の速度方程式は一般に：

$$
\frac{d\mathbf{x}}{dt} = \mathbf{f}(\mathbf{x})
$$

となる．
<!-- BZ反応が周期的に変動するのは反応速度の中に **強い非線形項** が存在するため． -->

代表的なモデル
- オレゴネータ(Oregonator)モデル
- Field–Körös–Noyes（FKN）機構
- 単純化すれば2次元の非線形常微分方程式でも振動が出る．

## Oregonatorモデル

- $x(t)$：中間体 HBrO$_2$ の濃度
- $y(t)$：還元剤の濃度
- $z(t)$：触媒の酸化状態に対応する変数

ここでは次の簡略版Oregonatorモデルを考える．

$$
\begin{cases}
\dot x = q y - x y + x(1 - x),
\\
\dot y = - q y - x y + f,
\end{cases}
$$

ここに
- $q>0$：反応比
- $f>0$：供給（定数）
- 非線形項は$xy$, $x^2$の2つ
- 適切なパラメタのもとでは，静止せずに閉じた軌道（リミットサイクル）に沿って動く．

## 平衡点の解析

平衡点は$\dot x=\dot y=0$で求まるが，ここでは表式が複雑になるため数値的に求める．

```python
import mpmath as mp

q, f = 0.002, 1.1
f1 = lambda X: q*X[1] - X[0]*X[1] + X[0]*(1-X[0])
f2 = lambda X: -q*X[1] - X[0]*X[1] + f

fpsol = mp.findroot([lambda X: f1(X), lambda X: f2(X)], (0.2, 1.0))
```

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

q = 0.002
f = 1.1

def oregonator_rhs(x, y):
    dx = q*y - x*y + x*(1-x)
    dy = -q*y - x*y + f
    return dx, dy
```


```python
def rk4_oregonator(f, x0, y0, t):
    X = np.empty_like(t); Y = np.empty_like(t)
    X[0], Y[0] = x0, y0
    h = t[1]-t[0]

    for k in range(len(t)-1):
        k1x, k1y = f(X[k], Y[k])
        k2x, k2y = f(X[k] + h*k1x/2, Y[k] + h*k1y/2)
        k3x, k3y = f(X[k] + h*k2x/2, Y[k] + h*k2y/2)
        k4x, k4y = f(X[k] + h*k3x,   Y[k] + h*k3y)
        X[k+1] = X[k] + (h/6)*(k1x + 2*k2x + 2*k3x + k4x)
        Y[k+1] = Y[k] + (h/6)*(k1y + 2*k2y + 2*k3y + k4y)

    return X, Y
```

計算の実施
```python
t = np.linspace(0, 200, 10000)
x0, y0 = 0.2, 1.2
X, Y = rk4_oregonator(oregonator_rhs, x0, y0, t)
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

![「研究内容」末松　J.　信彦](https://www.isc.meiji.ac.jp/~suematsu/research/pattern.html)

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

### 課題

1. Oregonator の2変数モデルを使って、異なる初期値から出発しても同じリミットサイクルに収束することを示せ。
2. パラメタ (q) と (f) を変化させ、安定なリミットサイクルが消滅・出現する様子を調べよ（ホップ分岐の観察）。
3. （発展）数値的に周期 (T) を推定せよ（ゼロクロス法など）。


