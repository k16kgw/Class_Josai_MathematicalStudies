第3回（ロジスティック成長）の体裁に合わせつつ，**第4回「変分原理と物理モデル」**の講義資料案を作りました．第3回で扱った「安定・不安定（ボウルの比喩）」を**ポテンシャル地形**に結び付け，そこから**最小作用の原理 → 単振り子**へ自然に橋渡しします．数式はそのままコピペで使えるようにしてあります．

---

# 変分原理と物理モデル

到達目標

- **最小作用の原理**の意味と，**オイラー＝ラグランジュ方程式**の導出を理解する．
- **単振り子モデル**を通じて**線形近似**（小振幅）と**非線形性の影響**（大振幅での周期変化）を理解する．

キーワード

- 最小作用の原理，オイラー＝ラグランジュ方程式
- ポテンシャル極値と安定性，位相平面
- 単振り子，線形化，非線形振動，周期の振幅依存

準備

1. anaconda を使用し，<span style="color:red">jupyter lab</span> を起動する．
2. `Documents（書類）/mathematical_studies` フォルダを開き，新規ノートブックを作成．
3. ファイル名を `4_{学籍番号}_{氏名}.ipynb` に変更（例：`4_SI25999_香川渓一郎.ipynb`）．

## 最小作用の原理

質点はその位置（配置）$x$に応じてポテンシャルエネルギー（$V(x)$）が定義される．

- ポテンシャルエネルギー：仕事を生み出す潜在的な「能力」
  - **極小**（$V^{\prime\prime}(x^*)>0$） $\Rightarrow$ **安定平衡**．
  - **極大**（$V^{\prime\prime}(x^*)<0$） $\Rightarrow$ **不安定平衡**．

- 現実として存在し得るのは「作用」が極小（最小）である状況であると理解できる．$\cdots$**最小作用の原理**

### 日常の最小作用

- **光の進み方：屈折・反射（フェルマーの原理＝最短時間）**
  - 2点間を進む光は，<span style="color:red">到達時間が最小</span>になる経路を進行する．
  - ガラスに斜めに入る光が屈折する（スネルの法則）．
- **石けん膜（最小面積）**
  - 針金枠に石けん膜を張ると，膜は<span style="color:red">面積が最小</span>の形に張る．
- **ぶら下がった鎖・電線の形（カテナリー曲線）**
  - 与えられた長さのもとで<span style="color:red">重力による位置エネルギーの総和が極小</span>になる形になる．

## スネルの法則

屈折率が $n_1, n_2$ の二媒質の境界に，光が入射角 $\theta_1$ で入ると，屈折角 $\theta_2$ は次で決まる．

$$
n_1 \sin\theta_1 = n_2 \sin\theta_2
$$

屈折率は$n := \dfrac{c}{v}$で定義される．
ここに$c$は真空中の光速，$v$は媒質中の位相速度．

遅い媒質へ入るときは法線側へ曲がり，速い媒質へは法線から離れる．

### 物理的な意味

- 光は「**到達時間が最小**」となる経路を選ぶ（**フェルマーの原理**）．
- 時間は次で得られる．

  $$
  T = \int \frac{ds}{v} = \frac{1}{c}\int n\ ds
  $$

### フェルマーの原理からの導出（最短時間）

上下で屈折率が一定の二層（上：$n_1$，下：$n_2$），境界は $y=0$．
出発点 ($x_A,y_A>0$) から到着点 ($x_B,y_B<0$) へ，境界上の通過点 $x$ に関して変分して**時間 $T(x)$** を最小化する．

- 距離と時間

  $$
  T(x)=\frac{1}{c}\Bigl[n_1 \sqrt{(x-x_A)^2+y_A^2} + n_2 \sqrt{(x_B-x)^2+y_B^2}\Bigr].
  $$

- 極値条件 ($\frac{dT}{dx}=0$) を計算すると

  $$
  &\quad
  \frac{n_1 (x-x_A)}{\sqrt{(x-x_A)^2+y_A^2}}
  =
  \frac{n_2 (x_B-x)}{\sqrt{(x_B-x)^2+y_B^2}}
  \\
  &\Leftrightarrow
  n_1\sin\theta_1 = n_2\sin\theta_2,
  $$

  となり，**スネルの法則**を得る．

<!-- 
### 波としての見方（ホイヘンス原理）

波面が境界に斜め入射すると，**境界上の位相連続**（波面の“縁”が繋がる）と**接線方向の波長連続**から，

$$
\frac{\sin\theta_1}{\lambda_1} = \frac{\sin\theta_2}{\lambda_2}
\quad\text{かつ}\quad
\frac{\lambda}{v}=\text{一定周波数}
;\Rightarrow;
n_1\sin\theta_1=n_2\sin\theta_2.
$$
-->
<!-- 
# 5) よく出る派生事項

* **全反射と臨界角**（(n_1>n_2) のとき）
  $$ n_1 \sin\theta_c = n_2 \quad\Rightarrow\quad \theta_c = \arcsin!\left(\frac{n_2}{n_1}\right). $$
  入射角 (\theta_1>\theta_c) で屈折光は出ず，**全反射**．
* **逆問題の対称性**
  行き来を反転しても同じ式．**可逆性**がある．
* **分散**（色によって曲がり方が違う）
  (n=n(\lambda))．プリズムで**青（短波長）ほど大きく曲がる**．

---

# 6) ベクトル表現（計算で便利）

入射単位ベクトル (\mathbf{s}_1)，法線単位ベクトル (\mathbf{n})（屈折面から入射側へ向き）を用いると，屈折方向 (\mathbf{s}_2) は

* 接線成分の連続：
  $$
  n_1,\mathbf{s}*{1,\parallel} = n_2,\mathbf{s}*{2,\parallel}
  $$
* 正規化 (|\mathbf{s}_2|=1)
  から一意に決定（3D レイトレーシング実装で使用）．
---

# 7) 例（数値感覚）

水面（空気 (n_1\approx 1.00) → 水 (n_2\approx 1.33)）．
入射角 (\theta_1=45^\circ) なら

$$
\sin\theta_2 = \frac{n_1}{n_2}\sin\theta_1 \approx \frac{1.00}{1.33}\times 0.707 \approx 0.532
;\Rightarrow;
\theta_2 \approx 32.1^\circ.
$$

**法線側に折れ曲がる**のが確認できます．
-->

### 応用例

- 眼鏡・レンズ設計，カメラ・望遠鏡
- 光ファイバー
- 異常震域
- 蜃気楼

### 屈折の経路の導出

準備
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'
```

```python
# 幾何（上が媒質1，下が媒質2，境界は y=0とする）
xA, yA = -2.0,  2.0   # 出発点 A (上側)
xB, yB =  3.0, -1.5   # 到着点 B (下側)

# 屈折率（例：空気→水）
n1 = 1.00   # 上側
n2 = 1.33   # 下側
c  = 3e8    # 真空中の光速（定数．相対比較なので消えてもOK）
```

```python
def travel_time(x):
    # 境界上の通過点 P=(x, 0)
    AP = np.hypot(x - xA, 0 - yA)  # 距離
    PB = np.hypot(xB - x, yB - 0)
    return (n1*AP + n2*PB) / c

# 点列を走査
xs = np.linspace(min(xA, xB) - 3, max(xA, xB) + 3, 500)
Ts = np.array([travel_time(x) for x in xs])

fig, ax = plt.subplots(figsize=(5,3))
ax.plot(xs, Ts)
ax.set_xlabel("境界上の位置 x")
ax.set_ylabel("到達時間 T(x) [任意単位]")
ax.set_title("到達時間 T(x) のプロファイル（凸な1谷型）")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

ここの$T$が最小になる$x$を見つける．
```python
i_min = np.argmin(Ts) # Tsの中で最小の値を取るindexを取得する
x_star = xs[i_min]
T_star = Ts[i_min]

AP = np.hypot(x_star - xA, yA)   # A側の経路の長さ
PB = np.hypot(xB - x_star, -yB)  # yB<0 を考慮した，B側の経路の長さ
sin1 = abs(x_star - xA) / AP
sin2 = abs(xB - x_star) / PB

lhs = n1 * sin1
rhs = n2 * sin2
```

経路を示す図の出力
```python
fig, ax = plt.subplots()
# 境界
ax.axhline(0, color='k', lw=1)
# 点と最適通過点
ax.plot([xA], [yA], 'o', label='A', color='blue')
ax.plot([xB], [yB], 'o', label='B', color='orange')
ax.plot([x_star], [0], 'o', label='P* (最適)', color='red')

# 経路
ax.plot([xA, x_star], [yA, 0], '-', color='blue')
ax.plot([x_star, xB], [0, yB], '-', color='orange')

ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()
```

## オイラー＝ラグランジュ方程式（Euler--Lagrange equation）

### 設定

- 一般化座標 $q$（座標の取り方に依らない）
- 一般化運動量 $\dot q$
- ポテンシャルエネルギー $V(q)$：$q$が大きいほど大きくなる
- 運動エネルギー $T(\dot q)$：$\dot q$が大きいほど大きくなる

- **ラグランジアン**

  $$
  L(q,\dot q,t)=T(\dot q)-V(q)
  $$
  
- **作用**

  $$
  S[q]=\int_{t_0}^{t_1} L\bigl(q(t),\dot q(t),t\bigr)dt
  $$

  経路 ($q(t)$) 全体に対し「<span style="color:red">どれだけ無駄なく動けたか</span>」を測る量．

### オイラー＝ラグランジュ方程式

ある時刻$t_0$で位置$q(t_0)=q_0$，$t_1$で$q(t_1)=q_1$を満たす何らかの現象を考える．
ここで**最小作用の原理**に基づき，作用を最小化する条件の下で得られるオイラー＝ラグランジュ方程式に従って現象が起きていると考える．
$q(t_0)=q_0$，$q(t_1)=q_1$を満たす現象$q(t)$が得られているとするとき，現象$q(t)$を僅かにずらしたときの作用の変化$\delta S$について$\delta S=0$を満たす$q$の方程式が求めるオイラー＝ラグランジュ方程式である．

任意の滑らかな関数$\eta(t)$で$\eta(t_0)=\eta(t_1)=0$を満たすものを取り，微小量$\varepsilon>0$について$q_\varepsilon(t) = q(t)+\varepsilon\eta(t)$とする．
このとき$q(t)$から$q_\varepsilon(t)$へずらしたときの$S[q]$の変化量$\delta S$（**変分**と呼ぶ）は

$$
\delta S &= S[q_\varepsilon] - S[q]
\\
&= \int_{t_0}^{t_1} L\bigl(q_\varepsilon(t),\dot q_\varepsilon(t),t\bigr)dt - \int_{t_0}^{t_1} L\bigl(q(t),\dot q(t),t\bigr)dt
$$

<!-- 
- **単位**：作用の単位はエネルギー×時間（J·s）．量子論で出てくるプランク定数 (h) も同じ単位．
 -->
<!-- 
## 0.5 よくある誤解と注意（3分）

- **「最小」じゃなくて** **「停留（極値）」**
  $$\delta S=0$$ は最小とは限らず，極大・鞍点の可能性もある（ただし典型的な力学では最小のことが多い）．

- **端点は固定**（今日の基本形）
  端点も自由なときは**自然境界条件**が別途出る（本編で一言触れる）．

- **一般化座標 (q)**
  直線座標だけでなく，角度や曲面上の座標でもよい（だから**拘束に強い**）．
-->

## 単振り子

### ラグランジアン

- 長さ $\ell$，質量 $m$，角度 $\theta$（下向き $0$）とする．
- 運動エネルギー・ポテンシャルエネルギー・ラグランジアンは次のようになる．

  $$
  T=\frac12 m\ell^2\dot\theta^2,\quad
  V=mg\ell(1-\cos\theta),\quad
  L=T-V.
  $$

### 運動方程式（E-L）

$$
\frac{d}{dt}(m\ell^2\dot\theta)+mg\ell\sin\theta=0
\quad\Rightarrow\quad
\boxed{\ \ddot\theta+\omega_0^2\sin\theta=0,\ \omega_0=\sqrt{g/\ell}\ }.
$$

### 平衡と安定・不安定

- ポテンシャル
  $$
  V(\theta)=mg\ell(1-\cos\theta)
  $$
  の極値：

  - (\theta=0)（極小）→ **安定平衡**（ボウルの底）．
  - (\theta=\pi)（極大）→ **不安定平衡**（伏せたボウル頂上）．

### エネルギー保存・位相平面

- 保存エネルギー

  $$
  E=\frac12 m\ell^2\dot\theta^2+mg\ell(1-\cos\theta).
  $$

- 位相平面 ($\theta,\dot\theta$)：

  - $E<2mg\ell$：**振動解**（閉曲線）．
  - $E=2mg\ell$：**セパラトリクス**（$\theta=\pi$ で停止／反転の境界）．
  - $E>2mg\ell$：**回転解**（一方向に回り続ける）．

### 線形近似（小振幅）

- $|\theta|\ll1$ で $\sin\theta\approx\theta$：

  $$
  \ddot\theta+\omega_0^2\theta=0,\quad
  \theta(t)=A\cos(\omega_0 t+\phi),\quad
  T_{\text{lin}}=\frac{2\pi}{\omega_0}=2\pi\sqrt{\frac{\ell}{g}}.
  $$

- **等時性**：周期は振幅に依らない（近似内）．

### 非線形性の影響（大振幅：周期は伸びる）

- 厳密周期は**完全楕円積分**で

  $$
  T(\theta_{\max}) = 4\sqrt{\frac{\ell}{g}};K(k),
  \quad
  k=\sin\frac{\theta_{\max}}{2},
  $$

  ここで (K) は第1種完全楕円積分．

- 小振幅展開：

  $$
  T(\theta_{\max}) \approx T_{\text{lin}}
  \left(1+\frac{1}{16}\theta_{\max}^2+\frac{11}{3072}\theta_{\max}^4+\cdots\right).
  $$

  ⇒ 振幅が大きいほど周期は**長く**なる．

### 線形モデルと非線形モデルの比較

準備
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'
```

パラメタと時間軸
```python
g = 9.8
ell = 1.0
w0 = np.sqrt(g/ell)

t_max = 12
dt = 1e-3
t = np.arange(0, t_max, dt)
```

小角近似の解
```python
def theta_lin(t, A, phi=0.0):
    return A*np.cos(w0*t + phi)
```

非線形の数値解（4次の古典RK）
```python
def pend_rhs(theta, omega):
    return omega, -(g/ell)*np.sin(theta)

def rk4(theta0, omega0):
    th = np.empty_like(t); om = np.empty_like(t)
    th[0], om[0] = theta0, omega0
    for k in range(len(t)-1):
        h = dt
        k1_th, k1_om = pend_rhs(th[k], om[k])
        k2_th, k2_om = pend_rhs(th[k]+0.5*h*k1_th, om[k]+0.5*h*k1_om)
        k3_th, k3_om = pend_rhs(th[k]+0.5*h*k2_th, om[k]+0.5*h*k2_om)
        k4_th, k4_om = pend_rhs(th[k]+h*k3_th, om[k]+h*k3_om)
        th[k+1] = th[k] + (h/6)*(k1_th + 2*k2_th + 2*k3_th + k4_th)
        om[k+1] = om[k] + (h/6)*(k1_om + 2*k2_om + 2*k3_om + k4_om)
    return th, om

A = np.deg2rad(30)  # 30°
theta_nl, omega_nl = rk4(A, 0.0)
theta_l = theta_lin(t, A, 0.0)
```

比較可視化
```python
fig, ax = plt.subplots(figsize=(5,3))
ax.plot(t, theta_l, label="小角近似", lw=1.5)
ax.plot(t, theta_nl, label="非線形（数値）", lw=1.5, alpha=0.9)
ax.set_xlabel("時間 [s]"); ax.set_ylabel(r"角度 $\theta(t)$ [rad]")
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout(); fig.savefig("./4_pend_compare.png", dpi=300); plt.show()
```

周期の推定（ゼロクロス法）
```python
def estimate_period(theta, t):
    # 0 を上向きに通過する時刻を拾う
    idx = np.where((theta[:-1] < 0) & (theta[1:] >= 0))[0]
    if len(idx) < 2: return np.nan
    t1 = t[idx[0]]; t2 = t[idx[1]]
    return t2 - t1

for deg in [5, 15, 30, 60]:
    A = np.deg2rad(deg)
    th, om = rk4(A, 0.0)
    T_est = estimate_period(th, t)
    print(f"初期振幅 {deg:>2}°: 推定周期 ≈ {T_est:.4f} s, 小角 {2*np.pi/w0:.4f} s")
```
<!-- 
期待される出力：振幅が大きいほど推定周期が**小角周期より長く**なる．
 -->

<!-- 
## 石鹸膜の形状

### 設定

- 2つの平行な円環を結ぶ石鹸膜の形状を求める．
- 半径 $R$ の円環ワイヤを $z=\pm \tfrac{L}{2}$ に平行に固定する．
- 両円環を結ぶ石鹸膜は**回転対称**とみなせるので，半径 $r=r(z)$ の回転面として表せる．
- 面積を$\mathcal{A} = \mathcal{A}[r]$は

$$
\mathcal{A}[r]= 2\pi \int_{-L/2}^{L/2} r(z),\sqrt{1+r'(z)^2}dz,
\qquad r!\left(\pm \tfrac{L}{2}\right)=R.
$$

### 解説（導出）

汎関数の被積分関数を

$$
\mathcal{L}(r,r')= r,\sqrt{1+r'^2}
$$

とおく．
$z$に陰にしか依存しないので，**ベルトラミの恒等式**が使える：

$$
\mathcal{L} - r'\frac{\partial \mathcal{L}}{\partial r'} = \text{const} = C.
$$

計算すると

$$
r\sqrt{1+r'^2} - r',\frac{r,r'}{\sqrt{1+r'^2}}
= \frac{r}{\sqrt{1+r'^2}} = C.
$$

よって

$$
\sqrt{1+r'^2}=\frac{r}{C}
\quad\Rightarrow\quad
r'^2 = \left(\frac{r}{C}\right)^2-1.
$$

変数分離して積分すると

$$
\int \frac{dr}{\sqrt{r^2-C^2}} = \int \frac{dz}{C}
\quad\Rightarrow\quad
\operatorname{arcosh}!\left(\frac{r}{C}\right) = \frac{z-z_0}{C}.
$$

したがって解は

$$
\boxed{\ r(z)= a,\cosh!\left(\frac{z-z_0}{a}\right)\ }, \qquad a(=C)>0.
$$

これが**カテノイド（catenoid）**：回転最小曲面の代表例．
境界条件 $r(\pm L/2)=R$ を用いると，（中心を対称に取れば $z_0=0$）

$$
\boxed{\ R = a\cosh\left(\frac{L}{2a}\right)\ } \quad (\star)
$$

を満たす $a$ が決まる．

### 面積の閉形式

$r=a\cosh(z/a)$, $r'=\sinh(z/a)$, $\sqrt{1+r'^2}=\cosh(z/a)$ より

$$
\mathcal{A}*{\text{cat}}
= 2\pi \int*{-L/2}^{L/2} a,\cosh^2!\left(\frac{z}{a}\right) dz
= 2\pi a^2\left[\sinh u \cosh u + u\right],
\quad u:=\frac{L}{2a}.
$$

境界式 ($\star$) で $R=a\cosh u$ を満たす $u>0$ が存在すれば，カテノイド解が存在．

### 重要な物理的含意（臨界間隔）

同半径 $R$ の2円環に対して，**カテノイドが存在するのは間隔が十分小さいときだけ**．
存在条件は上の境界式 $R=a\cosh u$（$a=L/(2u)$）の**最小可到達値**を考えればよい：

$$
R(u)=\frac{L}{2u}\cosh u \quad\text{の最小点で}\quad \frac{dR}{du}=0
\ \Rightarrow\ \boxed{\ \coth u = u\ }.
$$

この方程式の解は $u_\ast \approx 1.19968$．よって

$$
\boxed{\ \frac{L}{R}\le \frac{2u_\ast}{\cosh u_\ast} \approx 1.32549\ }
\quad\Bigl(\ \text{同値}\ \ \frac{L}{2R}\lesssim 0.66274\ \Bigr).
$$

この**臨界比**を超えると，カテノイドは存在しなくなり，石鹸膜は**2枚の円板**に「はじけて」移行します（実験で観察可能）． 
-->

<!-- 
## 演習問題

### 演習1（基礎）：「最小面積 ⇒ 最小曲面方程式」を導け

**問題**：グラフ面 (z=f(x,y)) の面積
$$
\mathcal{A}[f]=\iint_\Omega \sqrt{1+f_x^2+f_y^2},dx,dy
$$
の停留条件（E-L方程式）を導き，次を示せ：
$$
\boxed{\ \frac{\partial}{\partial x}!\left(\frac{f_x}{\sqrt{1+f_x^2+f_y^2}}\right)
+\frac{\partial}{\partial y}!\left(\frac{f_y}{\sqrt{1+f_x^2+f_y^2}}\right)=0\ }.
$$
**解答（略解）**：被積分関数 (\mathcal{L}=\sqrt{1+f_x^2+f_y^2}) に対して
(\partial \mathcal{L}/\partial f=0)，(\partial \mathcal{L}/\partial f_x=f_x/\sqrt{1+|\nabla f|^2}) 等を用い，2変数版E-Lを適用する．

### 演習2（中級）：「回転対称 ⇒ カテノイド」を自力で出す

**問題**：例題の回転面積汎関数
(\mathcal{A}[r]=2\pi\int r\sqrt{1+r'^2},dz) に対して，**ベルトラミの恒等式**を用いて
(\ r(z)=a\cosh!\bigl((z-z_0)/a\bigr)) を導け．
境界条件 (r(\pm L/2)=R) を課し，((\star)\ R=a\cosh(L/(2a))) を得よ．

**ヒント**：(\mathcal{L}) が (z) に陰にしか依存しない→
(\mathcal{L}-r'\partial \mathcal{L}/\partial r'= \text{const}) を使う．

### 演習3（上級）：「臨界間隔」を導け

**問題**：同半径 (R) の2円環を距離 (L) で向かい合わせたとき，
カテノイドが存在するための**臨界比** (L/R) を導け．
**解答（要点）**：
(R(u)=(L/2u)\cosh u) の最小を満たす (u) は (\coth u = u)．
その解 (u_\ast\approx 1.19968) を用いて
(\displaystyle (L/R)*{\max}= 2u*\ast/\cosh u_\ast \approx 1.32549).

### 演習4（比較）：「カテノイド vs 2枚の円板」のどちらが面積最小？

**問題**：カテノイドが存在する範囲で，
(\mathcal{A}*{\text{cat}}=2\pi a^2(\sinh u \cosh u + u)) と
(\mathcal{A}*{\text{disk}}=2\cdot \pi R^2) を比較し，
ある (L/R) 以上で **2枚の円板のほうが面積が小さくなる**ことを数値で確かめよ（転移点の概算を出せ）．
**ヒント**：(R=a\cosh u)，(L=2au)．パラメタ (u) で両者を評価し，差を可視化．

### 演習5（応用・数値）：「数値撮影で再現」

**問題**：任意の (R,L) を与えて，以下を計算・可視化せよ．

1. (u) を**数値的に**解く（(R=(L/2u)\cosh u)）．
2. (a=L/(2u)) を求め，(r(z)=a\cosh(z/a)) を (z\in[-L/2,L/2]) で描画．
3. 面積 (\mathcal{A}_{\text{cat}}) と (2\pi R^2) を比較．
   **ヒント**：数値解法は二分法やニュートン法で十分．図は (r(z)) を回転体として描くと直感的．

### 演習6（理論）：「平均曲率ゼロ」の確認

**問題**：最小曲面は**平均曲率 (H=0)** を満たす．
カテノイド (r(z)=a\cosh(z/a)) が（適切なパラメータ表示で）(H=0) を満たすことを示せ．
**ヒント**：回転曲面の平均曲率公式を用いるか，最小曲面方程式（演習1）を円筒座標に書き替えて確認．
 -->

## まとめ

| 観点        | 要点                                                                                                                  |
| --------- | --------------------------------------------------------------------------------------------------------------------- |
| 作用極値原理    | $\delta S=0\ \Rightarrow\ \frac{d}{dt}\Big(\frac{\partial L}{\partial \dot q}\Big)-\frac{\partial L}{\partial q}=0$ |
| 単振り子の力学   | $\ddot\theta+\frac{g}{\ell}\sin\theta=0,\ \ E=\tfrac12 m\ell^2\dot\theta^2+mg\ell(1-\cos\theta)$                |
| 安定・不安定    | $V''>0$で安定，$V''<0$で不安定
| 線形 vs 非線形 | 小角：等時性；大振幅：周期は $T(\theta_{\max})=4\sqrt{\ell/g},K(\sin(\theta_{\max}/2))$ で増大                           |
| 変分の利点     | 拘束・座標変換に強い．対称性→保存則の見通しが良い．                                                                         |

### 提出課題

1. 本日作成したipynbファイルをWebClassの「第4回課題」から提出せよ．

<!-- 
```{note}
**演習1（紙と鉛筆）**  
(1) 端点固定の仮定でオイラー＝ラグランジュ方程式を導出せよ．  
(2) ラグランジアン \(L=\tfrac12 m\dot x^2 - V(x)\) からニュートン方程式を導け．

**演習2（単振り子：線形化）**  
(1) \(\sin\theta\approx\theta\) で \(\ddot\theta+\omega_0^2\theta=0\) を得よ．  
(2) 小角周期 \(T_{\text{lin}}=2\pi\sqrt{\ell/g}\) を求めよ．

**演習3（非線形周期）**  
(1) エネルギー保存式から \(dt=\sqrt{\ell/(2g)}\,\frac{d\theta}{\sqrt{\cos\theta-\cos\theta_{\max}}}\) を導き，周期の積分形を示せ．  
(2) \(\theta_{\max}=10^\circ,30^\circ,60^\circ\) で，数値解から周期を見積もり，小角周期との差を議論せよ．

**演習4（発展：自然境界条件）**  
端点の位置ではなく「接線方向」が固定される最短曲線変分を考え，境界条件 \(\partial L/\partial \dot q=0\) がどのように現れるかを説明せよ．
```
 -->
