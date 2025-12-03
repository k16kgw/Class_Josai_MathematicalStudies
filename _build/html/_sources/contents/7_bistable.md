# 2種競争モデルと双安定性

到達目標
- 競争系における双安定現象を理解する．
- 相平面解析を行い安定性を評価する．

## 2種競争モデル

### 生態系の競争関係

- 一つの資源を2種以上の生物同士が奪い合う状況を扱う．
- 例：
  - 2種類の細菌が同じ栄養源を奪い合う．
  - 草原で2種の植物が日光を取り合う．
  - 2つの企業が同じ市場でシェア（客）を奪い合う．

→ このような**競争系**では，**一方が優勢になると他方が衰退する**現象が生じる．

### 競争排除則（Gauseの法則）

- 同じニッチ（生態的地位，まとまった環境要因）を占める2種は，長期的に共存できない．
- ただし，パラメータによっては「共存」や「多安定」もあり得る．

→ 代表例として**双安定現象**（bistability）を扱う．

### モデルの定式化

ロジスティック成長を基にモデルを考える．

単一種のロジスティック方程式は次のように導入される．

$$
\dfrac{dx}{dt} = r_1 x \left(1-\dfrac{x}{K_1}\right)
$$

これを2種 (x,y) が競合する場合に拡張する．

$$
\begin{cases}
    \dfrac{dx}{dt} = r_1 x \left(1-\dfrac{x}{K_1}-\dfrac{\alpha y}{K_1}\right)
    \\
    \dfrac{dy}{dt} = r_2 y \left(1-\dfrac{y}{K_2}-\dfrac{\beta x}{K_2}\right)
\end{cases}
$$

ここに
- $r_i$：種 $i$ の成長率
- $K_i$：種 $i$ の環境収容量
- $\alpha>0$：種1が種2に及ぼす競争効果（資源利用の強さ）
- $\beta>0$：種2が種1に及ぼす競争効果（資源利用の強さ）

**解釈**
- $x$ が増えると $y$ の成長が抑えられる（$\beta x$）
- $y$ が増えると $x$ の成長が抑えられる（$\alpha y$）
- 強い競争ではどちらかが勝ち残り，弱い競争では共存も可能．

---

## 平衡点とその分類

### 平衡条件

$\dfrac{dx}{dt}=\dfrac{dy}{dt}=0$より

$$
\begin{cases}
    x=0 \text{ または } 1-\dfrac{x}{K_1}-\dfrac{\alpha y}{K_1}=0,
    \\
    y=0 \text{ または } 1-\dfrac{y}{K_2}-\dfrac{\beta x}{K_2}=0.
\end{cases}
$$

ここで簡単のため $r_1=r_2=1,\ K_1=K_2=1$ として考えると：

$$
\begin{cases}
    \dfrac{dx}{dt} = x(1 - x - \alpha y)
    \\
    \dfrac{dy}{dt} = y(1 - y - \beta x)
\end{cases}
$$

ここで第1・2式の右辺を$f_1, f_2$と置く．つまり

$$
&f_1 = x(1 - x - \alpha y)
\\
&f_2 = y(1 - y - \beta x)
$$

### 平衡点の候補

1. $E_1=(0,0)$：両方絶滅
2. $E_2=(1,0)$：種1のみ生存
3. $E_3=(0,1)$：種2のみ生存
4. $E_4=\left(\dfrac{1-\alpha}{1-\alpha\beta}, \dfrac{1-\beta}{1-\alpha\beta}\right)$：共存状態（ただし存在条件あり）

### 共存条件

共存点が正（第1象限内）に存在するための条件は

$$
1-\alpha>0,\quad 1-\beta>0 \quad\Rightarrow\quad \alpha<1,\ \beta<1.
$$

→ 相互作用が弱ければ共存可能．

逆に $\alpha>1,\ \beta>1$ では，共存点は第1象限外になり，**排除**（双安定）が起きる．

## 相平面解析

### ヌルクライン

$\dfrac{dx}{dt}=0, \dfrac{dy}{dt}=0$を満たす$(x,y)$の関係式（曲線）を$xy$平面に描画することで，相平面内の軌道の大まかな挙動を予測することができる．
この曲線は

$$
&f_1 = x(1 - x - \alpha y) = 0
\\
&f_2 = y(1 - y - \beta x) = 0
$$

より

$$
&x = 0
\\
&y = - \dfrac{1}{\alpha} x + \dfrac{1}{\alpha}
\\
&y = 0
\\
&y = - \beta x + 1
$$

と表せる．この曲線を境に$\dfrac{dx}{dt}$と$\dfrac{dy}{dt}$の正負の値が変わり，$x,y$の増減が決まる．



### 数値的な解析

**設定**
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'
```

```{note}
**演習1**

次のコード内で2種競争モデルの右辺を計算する関数`competition_rhs`を実装せよ．
```python
alpha, beta = 1.3, 1.4  # 強い競争
def competition_rhs(x, y):
    """
    ここを埋める．
    """
    return dx, dydt
```
<!-- 
```python
alpha, beta = 1.3, 1.4  # 強い競争
def competition_rhs(x, y):
    dx = x*(1 - x - alpha*y)
    dy = y*(1 - y - beta*x)
    return dx, dy
```
 -->

**ベクトル場の描画**
```python
x = np.linspace(0, 1.5, 25)
y = np.linspace(0, 1.5, 25)
X, Y = np.meshgrid(x, y)
U, V = competition_rhs(X, Y)

fig, ax = plt.subplots(figsize=(5,4))
ax.quiver(X, Y, U, V, color='gray', angles='xy')
ax.set_xlabel("種1の個体数 x")
ax.set_ylabel("種2の個体数 y")
ax.set_title("2種競争モデルの方向場")
ax.grid(True, alpha=0.3)
plt.show()
```

**時間発展と軌道（オイラー法）**
```python
def euler_comp(f, x0, y0, t):
    X = np.empty_like(t); Y = np.empty_like(t)
    X[0], Y[0] = x0, y0
    h = t[1]-t[0]
    for k in range(len(t)-1):
        dxdt, dydt = f(X[k], Y[k])
        X[k+1] = X[k] + h*dxdt
        Y[k+1] = Y[k] + h*dydt
    return X, Y

t = np.linspace(0, 50, 2000)
x0, y0 = 0.2, 1.0
X, Y = euler_comp(competition_rhs, x0, y0, t)
```

**相平面上の軌道**
```python
fig, ax = plt.subplots(figsize=(5,4))
ax.plot(X, Y, 'b-', label='軌道')
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_xlabel("x (種1)")
ax.set_ylabel("y (種2)")
ax.set_title("競争系の相平面軌道（双安定）")
ax.legend(); ax.grid(True, alpha=0.3)
plt.show()
```

```{note}
**演習2**

いくつか異なる初期値について軌道を観察することで，$E_2=(1,0)$に収束する場合と$E_3=(0,1)$に収束する場合との条件の仮説を立てよ．
<!-- 
- 種1の初期値が大きいと $x\to 1, y\to 0$．
- 種2が優勢なら $x\to 0, y\to 1$．
- 中間値からの出発では，**どちらに偏るかが初期条件に依存**する．

→ これを**双安定**（bistability）と呼ぶ．
 -->
```

## 双安定現象

- 系が **2つの安定平衡点** を持つとき，そのどちらにも収束し得る現象をいう．
- 初期条件のわずかな違いで最終状態が異なる．

### 相図による視覚的理解

- 境界（セパラトリクス）を挟んで，軌道が異なる安定平衡点に吸い込まれる．
- 実世界では：
  - 競争種の**どちらが勝つか**
  - 化学反応での**スイッチング挙動**
  - 生物細胞の**分化** などに対応する．

### 平衡点の安定性の解析

$$
J=
\begin{pmatrix}
    1-2x-\alpha y   & -\alpha x
    \\
    -\beta y        & 1-2y-\beta x
\end{pmatrix}
$$

各平衡点に代入して固有値の符号を確認すると

| 平衡点 | 固有値の符号 | 安定性         | 意味            |
| :---: | :---:      | :-----------  | :------------- |
| (1,0) | 負・正      | 安定（種1のみ） | 種1が優勢        |
| (0,1) | 正・負      | 安定（種2のみ） | 種2が優勢        |
| (0,0) | 正・正      | 不安定        | 両者増える        |
| 共存点 | なし        | なし          | 強い競争のため消失 |

2つの安定平衡$(1,0)$と$(0,1)$が存在：**双安定系**．

```{note}
**演習3**

1. 弱い競争（$\alpha=0.5, \beta=0.6$）に変更し，共存点が出現することを確認せよ．
2. 初期条件を変えて軌道がどちらの平衡点へ向かうかを可視化せよ．
<!-- 3. （発展）ヤコビ行列の固有値を数値的に計算し，安定性を評価せよ． -->
```

---

## まとめ

| 項目    | 内容                               |
| ----- | -------------------------------- |
| モデル   | 2種競争モデル（ロジスティック型）                |
| 双安定性  | 初期条件により異なる安定平衡に収束                |
| 安定性解析 | ヤコビ行列の固有値の符号で判断                  |
| 相平面解析 | 境界を挟んで異なる方向に流れるベクトル場             |
| 応用例   | 種間競争，化学反応のスイッチング，神経モデル，細胞分化モデルなど |

