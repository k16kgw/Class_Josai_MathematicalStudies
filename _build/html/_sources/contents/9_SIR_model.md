# 感染症モデル（SIRモデル）

到達目標
- SIRモデルの各変数 $S,I,R$ が**現実の何を表すか**を明確に説明できる．
- SIRモデルの基本構造（不変集合・閾値・最終感染規模）を解析的に導出できる．
- 実データとモデルの比較を「観測モデル＋パラメータ同定（逆問題）」として定式化できる．
- モデルが実データを再現できない場合に，その原因を数理的に分類できる．

キーワード
- コンパートメントモデル（区分モデル）
- 正不変性・保存則・単調性
- 基本再生産数 $R_0$
- 最終感染規模方程式（final size equation）
- 観測モデル・報告率・逆問題（同定）
- モデル誤差（model discrepancy）

準備
1. anacondaを使用し，<span style="color:red">jupyter lab</span>を起動する．
2. `Documents（書類）/mathematical_studies`フォルダをダブルクリックで開き`+`をクリックして新しいファイルを作成する．
3. ファイル名を`9_{学籍番号}_{氏名}.ipynb`に変更する．例：`9_SI25999_香川渓一郎.ipynb`

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'
```

---

## 感染症を数理モデルで捉える

新型コロナウイルスの感染拡大をデータで見返す．

![感染者数](/contents/figs/9/count.png)

※ 日本経済新聞「チャートで見る日本の感染状況 新型コロナウイルス」2025年1月8日閲覧．[https://vdata.nikkei.com/newsgraphics/coronavirus-japan-chart/](https://vdata.nikkei.com/newsgraphics/coronavirus-japan-chart/)

感染症の拡大は個体間の接触によって起こる感染が連鎖することによって進む．

- 感染性をもつAさんと感受性をもつBさんが接触する．
  - 感染性：他人にウイルスを感染させる性質
  - 感受性：他人からウイルスを感染させられる性質（まだ感染しておらず，免疫を持たない）
- 一定の確率でBさんが感染する．

この現象を数理モデルで捉える．
- これは本質的に**確率過程**であり，厳密に扱うなら個体をノードとしたネットワーク上の感染伝播を考えることになる．
- ただしそのままでは次の困難が生じる．
  - 個体数が多すぎる．
  - 接触ネットワークが未知である．
  - データから推定する自由度が大きすぎる．
- そこで，「個体の詳細は捨て，集団全体を平均化して記述する」方向で数理モデル化を試みる．

感染症の拡大において，人々を次の3種類に大別する．
- 感受性を持つ人（Susceptible）
- 感染性を持つ人（Infectious）
- 除去された人（Removed）

これらの頭文字を取って<span style="color:red">SIRモデル</span>と呼ばれる．

### $S(t)$：感受性を持つ人の数

- まだ感染しておらず，**感染する可能性がある人**の数．
- “免疫がない（十分でない）人” と読み替えてよい．
- 現実には，ワクチン・既感染・年齢層などで「感受性」は層構造を持つが，SIRでは**全部まとめて1つ**にする．

### $I(t)$：感染性を持つ人の数

- 他者へ**ウイルスを移し得る状態にいる人**の数．
- 検査で陽性になった人数ではないことに注意．
- 「無症状」でも感染性があればこれに分類される．
- 観測されるデータ（報告感染者数）は多くの場合，感染性を持つ人の一部に過ぎない．

### $R(t)$：除去された人（Removed）

- 以下のように**感染連鎖から外れた人**をまとめて表すことが多い．
  - 回復して免疫を得た（再感染）
  - 死亡した
  - 隔離された
  - 入院して接触しない

### 観測できる量

現実のデータは次のようなものがある：
- 日次新規感染者（incidence）
- 検査陽性者数（reported cases）
- 入院者数・重症者数
- 死亡者数

しかし，これらは $S,I,R$ のどれかを“直接”見ているわけではない．

---

## SIRモデルの導出

### モデルの仮定（現象論）

- 時刻 $t$ における感染機会の総数は，感受性者$S$と感染者$I$の**組み合わせ**に比例すると考える．
- 1人の感染者が単位時間に作る接触の回数×感染確率をまとめて $\beta$ と置く．
- 除去は感染者に比例し，比例定数を $\gamma$ と置く．

### SIRモデル

$$
\begin{aligned}
    \frac{dS}{dt} &= -\beta SI,
    \\
    \frac{dI}{dt} &= \beta SI-\gamma I,
    \\
    \frac{dR}{dt} &= \gamma I.
\end{aligned}
$$

- $\beta$：感染が起きる速さ（社会的接触・対策の影響が強い）
- $\gamma$：除去の速さ（回復・隔離の速さ）

---

## 力学系としての基本構造（解析の入口）

### 正不変性（非負性）

初期値が $S,I,R\ge 0$ なら，時間発展で負にならない．

→ 「人数を表しているのに負になる」という破綻は起こらない．

```{note}
初期値が

$$
S(0)\ge 0,\quad I(0)\ge 0,\quad R(0)\ge 0
$$

を満たすとき，解が存在する限り

$$
S(t)\ge 0,\quad I(t)\ge 0,\quad R(t)\ge 0\qquad (\forall t\ge 0)
$$

が成り立つ。すなわち

$$
\mathbb{R}^3_{\ge 0}=\{(S,I,R)\mid S\ge0,\ I\ge0,\ R\ge0\}
$$

は **正不変集合**である。
```

*証明*

ベクトル場

$$
F(S,I,R)=(-\beta SI,\ \beta SI-\gamma I,\ \gamma I)
$$

を考える．

**(i) 平面 $S=0$ 上**

$S=0$ なら

$$
\frac{dS}{dt}=-\beta\cdot 0\cdot I = 0.
$$

よって，軌道は$S=0$を横切って$S<0$側へ出ていくことはない（少なくとも$S$成分は負方向に進めない）．

```{note}
<span style="color:red">**演習1**</span>

平面 $I=0$，$R=0$ についても軌道が$I<0$，$R<0$へと出ていかないことを確かめよ．
```

<!-- 
**(ii) 平面 $I=0$ 上**

$I=0$ なら

$$
\frac{dI}{dt}=\beta S\cdot 0-\gamma\cdot 0 = 0.
$$

よって，軌道は $I=0$ を横切って $I<0$ 側へ出ていくことはない．

**(iii) 平面 $R=0$ 上**

$R$ の式は $R$ 自体を含まず

$$
\frac{dR}{dt}=\gamma I\ge 0 \quad (I\ge 0)
$$

だから，非負領域内では $R$ は減らず，$R<0$ に落ちることはない．

以上より，非負直交体の境界から外（負の方向）へ向かう成分が存在しないため，初期値が非負なら解は常に非負に留まる．□
 -->

### 保存則

$$
\frac{dS}{dt} + \frac{dI}{dt} + \frac{dR}{dt} = 0.
$$

より

$$
S(t)+I(t)+R(t)=N=\mathrm{const.}
$$

3変数から2変数に減らすことができる．

### 単調性

```{note}
次が成り立つ．

$$
\frac{dS}{dt} \le 0,
\\
\frac{dR}{dt} \ge 0.
$$

すなわち
- $S(t)$は単調減少
- $R(t)$は単調増加
```

*証明*

**(i) $S(t)$ は単調減少する．**

$$
\frac{dS}{dt}=-\beta S I \le 0
\qquad(\beta>0,\ S\ge0,\ I\ge0).
$$

従って $S(t)$ は単調非減少する．

<!-- - さらに (S(t)) が一定になるのは，区間上で $SI\equiv 0$（つまり $S\equiv 0$ または $I\equiv 0$）のときに限られる。 -->

**(2) $R(t)$ は単調増加する**

$$
\frac{dR}{dt}=\gamma I \ge 0
\qquad(\gamma>0,\ I\ge0).
$$

従って $R(t)$ は単調増加する．

<!-- - $R(t)$ が一定なのは $I\equiv 0$ のときに限られる。 -->

```{warning}
一般に$I(t)$は単調ではない．
また，$I(t)$は極大値を取るとしたら1回のみで，その後は単調減少する．（単峰性）
```

$$
\frac{dI}{dt} = I(\beta S-\gamma),
$$

であるから $I(t)$ の増減は $S(t)$ の大きさで決まる：

- $S(t)>\gamma/\beta$ のとき：$\frac{dI}{dt}>0$（増加）
- $S(t)=\gamma/\beta$ のとき：$\frac{dI}{dt}=0$（極値候補）
- $S(t)<\gamma/\beta$ のとき：$\frac{dI}{dt}<0$（減少）

ここで $S(t)$ は単調減少なので，$S(0)>\gamma/\beta$であれば $I(t)$ は最初増加し，ある時刻でピークを迎え，その後減少する（単峰）という形になる．

もちろん$S(0)<\gamma/\beta$であれば$I(t)$はピークを迎えずに減少する．

---

## 閾値理論

現実問題として，初期値では人口のほとんどは感受性者であり，そこに僅かな感染者がいる状況を考える．
まだ感染は拡大していないので，除去された人はいないものとする．

$$
S(0) &\approx N,
\\
I(0) &\ll 1,
\\
R(0) &= 0.
$$

### 初期段階の近似

初期は $S(t) \approx N$ とみなせるので

$$
\frac{dI}{dt} \approx (\beta N - \gamma) I.
$$

従って，

$$
R_0:=\frac{\beta N}{\gamma}
$$

とおくと

- $R_0>1$：感染拡大
- $R_0<1$：感染終息

を表す．

この$R_0$を<span style="color:red">基本再生産数</span>と呼ぶ．

- 接触が多い（大きい $\beta$）ほど$R_0$は大きい
- 隔離が早い（大きい $\gamma$）ほど$R_0$は小さい

→ 政策介入によって感染拡大か終息かを制御できる可能性を示唆．

```{note}
<span style="color:red">**演習2**</span>

感染を拡大させないためには$R_0$を小さくすることが肝要である．
$R_0$を小さくするためにできることには何があるか，5つ以上考えよ．
```

---

<!-- ## ピーク条件と最終感染規模方程式 -->
## 最終感染規模方程式（Final size equation）

<!-- ### 感染ピークの条件

ある時刻$t_{\text{peak}}$で感染者数$I$がピークを取るとき，$\frac{dI}{dt}(t_{\text{peak}})=0$であるから，
$$
\frac{dI}{dt}(t_{\text{peak}}) = I(t_{\text{peak}})\left(\beta S(t_{\text{peak}}) - \gamma\right)=0,
$$
より
$$
S(t_{\text{peak}})=\frac{\gamma}{\beta}=\frac{N}{R_0}.
$$ -->

<!-- ### 最終感染規模方程式（Final size equation） -->

$t\to \infty$で感染が収束（$I(t)\to 0$）したときの感受性者数$S(t)\to S_{\infty}$と除去された人の数$R(t)\to R_{\infty}$を調べる．

$\frac{dS}{dt}$と$\frac{dR}{dt}$の微分方程式より次の微分方程式が得られる．

$$
\frac{dS}{dR}=-\frac{\beta}{\gamma} S
= - R_0 \frac{S}{N}.
$$

この微分方程式は変数分離形式なので解ける．
保存則$S(t)+I(t)+R(t)=N$と$t\to\infty$の極限を考えれば，未知数を$S_{\infty}$のみにすることができ，次の<span style="color:red">最終感染規模方程式</span>を得る．

$$
\log S_\infty + \frac{\beta}{\gamma}(N-S_\infty) = \log S(0) + \frac{\beta}{\gamma}R(0).
\tag{FSE}
$$

ここで$R(0)=0$なので

$$
\frac{S_\infty}{S_0}
= \exp\left(-\frac{\beta}{\gamma}(N-S_\infty)\right),
$$

とできる．
これを解くことで$S_\infty$を求めることができるが，超越方程式なので数値的に解くことを試みる．

SIRモデルの時間発展を全て追わなくても，

- 最終的に感受性者がどれだけ残るか（$S_\infty$）
- 累積感染者（除去者）数がどれくらいか（$R_\infty=N-S_\infty$）

がこの方程式から求まる．

特に

- $\beta/\gamma = \frac{R_0}{N}$ が大きいほど $S_\infty$ は小さくなり，
- 累積感染者数 $R_\infty$ が大きくなる．

### 数値的に最終規模方程式を解く

```python
def fse_residual(S, S0, N, beta, gamma, R_init=0.0):
    """FSEの左辺-右辺（=0を解く）"""
    a = beta / gamma
    return np.log(S) + a*(N - S) - (np.log(S0) + a*R_init)

def solve_fse_bisection(S0, N, beta, gamma, R_init=0.0, tol=1e-10, max_iter=200):
    """
    二分法で S_inf を解く。
    S_inf は (0, S0] にある（流行でSは減るので）。
    """
    # loからhiの間で解を探す
    # 下端は log のため0にできないので小さい正数
    lo = 1e-15 * N
    hi = S0

    f_lo = fse_residual(lo, S0, N, beta, gamma, R_init)
    f_hi = fse_residual(hi, S0, N, beta, gamma, R_init)

    # 通常パラメータでは f(lo)>0, f(hi)<=0 となり根を挟む
    # もし挟めない場合は情報を出して止める
    if f_lo * f_hi > 0:
        raise ValueError(
            "二分法の区間が根を挟んでいません。"
            f" f(lo)={f_lo:.3e}, f(hi)={f_hi:.3e}. "
            "N,beta,gamma,S0 の設定を確認してください。"
        )

    for _ in range(max_iter):
        mid = 0.5*(lo+hi)
        f_mid = fse_residual(mid, S0, N, beta, gamma, R_init)

        if abs(f_mid) < tol or (hi-lo) < tol*max(1.0, hi):
            return mid

        # 根を挟む側を更新
        if f_lo * f_mid > 0:
            lo, f_lo = mid, f_mid
        else:
            hi, f_hi = mid, f_mid

    return 0.5*(lo+hi)
```

```python
N = 1
S0 = N - 0.0001   # 初期感染がわずか

beta = 0.30
gamma = 0.10
R_init = 0.0

S_inf = solve_fse_bisection(S0, N, beta, gamma, R_init=R_init)
R_inf = N - S_inf

print(f"S_inf ≈ {S_inf:.6f}")
print(f"R_inf ≈ {R_inf:.6f}")
print(f"最終感染割合 R_inf/N ≈ {R_inf/N:.6f}")
print(f"beta/gamma = {beta/gamma:.3f}")
```

最終規模を可視化
```python
N = 1.0       # 正規化して割合で見たい場合（N=1が便利）
S0 = N - 0.00001      # 初期ほぼ全員感受性（割合）
R_init = 0.0

R0_list = np.linspace(0.5, 5.0, 46)  # 0.5〜5.0
Sinf_list = []
Rinf_list = []

for R0 in R0_list:
    beta = R0  # gamma=1 とすれば beta/gamma=R0
    gamma = 1.0
    S_inf = solve_fse_bisection(S0, N, beta, gamma, R_init=R_init)
    Sinf_list.append(S_inf)
    Rinf_list.append(N - S_inf)

Sinf_list = np.array(Sinf_list)
Rinf_list = np.array(Rinf_list)

plt.figure(figsize=(7,3))
plt.plot(R0_list, Sinf_list, label=r"$S_\infty/N$")
plt.plot(R0_list, Rinf_list, label=r"$R_\infty/N = 1-S_\infty/N$")
plt.axvline(1.0, color="gray", lw=1, alpha=0.6)
plt.xlabel(r"$R_0=\beta/\gamma$")
plt.ylabel("ratio")
plt.title("基本再生産数に対する最終感染規模方程式の解")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
```

---

## SIRモデルの限界

実際の感染の動向では感染者数のピークは何度も現れる．

![感染者数](/contents/figs/9/count.png)

SIRモデルでの変数$I(t)$ある時刻$t$における感染者数の**総数**を記述しており，その日に報告された感染者数とは厳密には一致しない．
感染者がほぼ全て入院するものと仮定すれば，その日にいる患者数が最も$I(t)$の意味に近いだろう．

![患者数](/contents/figs/9/hospital.png)

※ 日本経済新聞「チャートで見る日本の感染状況 新型コロナウイルス」1月8日閲覧．[https://vdata.nikkei.com/newsgraphics/coronavirus-japan-chart/](https://vdata.nikkei.com/newsgraphics/coronavirus-japan-chart/)

患者数の動向を観ると何度かピークが現れていることが分かる．
しかし，SIRモデルでは$I(t)$は単峰性を持つことから，ピークは1度しか訪れない．

→ 現実のデータに合わせるにはモデルを改善する必要がある．
（**パラメタの値を変えるだけでは現実のデータには合わせられない**）

```{note}
<span style="color:red">**演習3**</span>

より現実に即したモデルにするには何を変えれば良いか．
5つ以上アイデアを挙げよ．
```
<!-- 
### モデルの拡張

- $\beta=\beta(t)$（時変：政策・行動変容・季節性）
- SEIR（潜伏期）
- 免疫減衰 SIRS
- 空間・ネットワーク
- 複数株
 -->

<!-- 
## 実データとの比較：観測モデル → 逆問題（同定）へ

観測は $I$ そのものではなく投影であり，だから逆問題になる．

### 観測モデルの例

例A：有病者（感染性者）に比例

$$
Y(t_k)=\rho I(t_k)+\eta_k.
$$

- $\rho$：報告率（未観測）
- $\eta_k$：ノイズ

例B：日次新規感染者（incidence）

SIRの感染発生率は

$$
\text{inc}(t)=\beta\frac{S(t)I(t)}{N}.
$$

日次データなら

$$
Y_k \approx \rho \int_{t_k}^{t_{k+1}}\beta\frac{S(t)I(t)}{N},dt + \eta_k.
$$

### 同定問題（最小二乗）

$$
\min_{\beta,\gamma,\rho,I_0}
\sum_k\left(Y_k-\hat Y_k(\beta,\gamma,\rho,I_0)\right)^2.
$$

ここで重要な注意：

- 初期指数増加だけでは $\beta-\gamma$ しか分かりにくい
- $\rho$ と $I$ は分離不能になりやすい
  → **同定不能性**が起きる

この“数学的限界”を意識した上で推定をする．

---

## 実装（実習）：データ読み込み → フィット → 残差

### データ読み込みと前処理

```python
df = pd.read_csv("data_daily_cases.csv")
df["date"] = pd.to_datetime(df["date"])
df = df.sort_values("date")

y = df["new_cases"].to_numpy()
t = np.arange(len(y), dtype=float)

# 7日移動平均で平滑化
w = 7
y_ma = np.convolve(y, np.ones(w)/w, mode="same")

plt.figure(figsize=(7,3))
plt.plot(df["date"], y, alpha=0.25, label="raw")
plt.plot(df["date"], y_ma, label="7-day MA")
plt.legend(); plt.grid(True, alpha=0.3)
plt.show()
``` -->

### SIRの数値解

```python
def sir_rhs(S, I, R, beta, gamma, N):
    dS = -beta*S*I/N
    dI = beta*S*I/N - gamma*I
    dR = gamma*I
    return dS, dI, dR

def euler_sir(S0, I0, R0, beta, gamma, N, t):
    S = np.empty_like(t, dtype=float)
    I = np.empty_like(t, dtype=float)
    R = np.empty_like(t, dtype=float)
    S[0], I[0], R[0] = S0, I0, R0

    h = t[1] - t[0]
    for k in range(len(t)-1):
        dS, dI, dR = sir_rhs(S[k], I[k], R[k], beta, gamma, N)
        S[k+1] = S[k] + h*dS
        I[k+1] = I[k] + h*dI
        R[k+1] = R[k] + h*dR

        # 数値誤差で負に落ちるのを避ける（教育上の安全策）
        S[k+1] = max(S[k+1], 0.0)
        I[k+1] = max(I[k+1], 0.0)
        R[k+1] = max(R[k+1], 0.0)

    return S, I, R
```

```python
N = 1_000_000
beta = 0.30
gamma = 0.10

S0 = N - 10
I0 = 10
R0 = 0

# 時間軸（例：0〜160日）
t = np.linspace(0, 160, 1601)  # 刻み幅 h=0.1

# ===== 数値計算 =====
S, I, R = euler_sir(S0, I0, R0, beta, gamma, N, t)

# （任意）日次新規感染者（incidence）
inc = beta * S * I / N

# ===== 可視化：S, I, R の時系列 =====
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(t, S, label="S(t) 感受性者")
ax.plot(t, I, label="I(t) 感染性者")
ax.plot(t, R, label="R(t) 除去者")
ax.set_xlabel("時刻 t")
ax.set_ylabel("人数")
ax.set_title(f"SIRモデル（Euler）  R0=beta N/gamma={beta*N/gamma:.2f}")
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
plt.show()

# ===== 可視化：日次新規感染者（incidence） =====
fig, ax = plt.subplots(figsize=(7, 3.0))
ax.plot(t, inc, label=r"inc(t) = $\beta S I / N$")
ax.set_xlabel("時刻 t")
ax.set_ylabel("新規感染者（モデル）")
ax.set_title("日次新規感染者（incidence）の時系列")
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
plt.show()
```


<!-- 
### 観測：incidence を当てる

```python
from scipy.optimize import least_squares

N = 1_000_000  # 有効人口（まず固定）

def model_incidence_euler(params):
    beta, gamma, rho, I0 = params
    S0, R0_ = N - I0, 0.0

    # データと同じ時刻グリッドで解く
    t = t_data.copy()
    S, I, R = euler_sir(S0, I0, R0_, beta, gamma, N, t)
    inc = beta*S*I/N
    return rho*inc

def residuals(params):
    return model_incidence_euler(params) - y_ma

x0 = np.array([0.3, 0.1, 0.2, 10.0])
lb = np.array([0.0, 0.0, 0.0, 1.0])
ub = np.array([5.0, 5.0, 10.0, 1e6])

res = least_squares(residuals, x0, bounds=(lb, ub))
beta_hat, gamma_hat, rho_hat, I0_hat = res.x

print("beta =", beta_hat)
print("gamma =", gamma_hat)
print("rho =", rho_hat)
print("I0 =", I0_hat)
print("R0 =", beta_hat/gamma_hat)
```

### フィットと残差の評価

```python
y_fit = model_incidence(res.x)

plt.figure(figsize=(7,3))
plt.plot(y_ma, label="data (MA)")
plt.plot(y_fit, label="SIR incidence fit")
plt.legend(); plt.grid(True, alpha=0.3)
plt.show()

plt.figure(figsize=(7,3))
plt.plot(y_ma - y_fit)
plt.title("residual (data - model)")
plt.grid(True, alpha=0.3)
plt.show()
```

ここで残差に系統的構造（複数波など）が見えれば，それは推定の失敗ではなく **モデル構造の不足**．

---

## まとめ

- §1で「$S,I,R$は何を数えているか」を固定した
  → 実データが直接 $I$ ではないことが分かった
- §2でモデルを導出し，パラメータ $\beta,\gamma$ の意味が明確になった
- §3〜§5で，解の全体像（閾値・ピーク・最終感染規模）を解析的に得た
- §6で「SIRは単峰になりやすい」構造的理由を示した
- §7〜§8で，観測モデルを入れて逆問題としてフィットし，残差からモデル誤差を評価した
- 次回：モデル構造の不足を埋める（$\beta(t)$，SEIR，SIRS，空間など） -->

---
<!-- 
## 演習（レベル設定：応用数学）

**演習1（解析）**
SIRから最終感染規模方程式(FSE)を導出し，与えられたR0についてS(∞)を数値的に解け．

**演習2（理論）**
S(t)の単調性を用いて，I(t)が高々1回しか極大を持たないことを示せ．

**演習3（逆問題）**
同定で得られた(beta, gamma)が初期値x0に依存するか調べ，局所最小の可能性を議論せよ．

**演習4（拡張）**
β(t)=β0(1+a sin(2πt/T))を導入し，複数波が出ることを確認せよ（数値実験）．
 -->

## 現実のデータを見る

[https://idsc.tmiph.metro.tokyo.lg.jp/diseases/flu/flu/](https://idsc.tmiph.metro.tokyo.lg.jp/diseases/flu/flu/)

```python
# 取得したデータ
y = np.array([
    161, 278, 416, 818, 1385, 1990, 2335, 4333,
    9926, 12133, 18707, 21608, 14947, 10495, 7440, 6739
], dtype=float)
t = np.arange(len(y), dtype=float)  # 週インデックス（0,1,2,...）

plt.figure(figsize=(7,3))
plt.plot(t, y, "o-", label="data (weekly total)")
plt.xlabel("week index")
plt.ylabel("cases")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
```

```python
def sir_rhs_frac(s, i, r, beta, gamma):
    ds = -beta*s*i
    di = beta*s*i - gamma*i
    dr = gamma*i
    return ds, di, dr

def euler_sir_frac(s0, i0, r0, beta, gamma, t):
    s = np.empty_like(t, dtype=float)
    i = np.empty_like(t, dtype=float)
    r = np.empty_like(t, dtype=float)
    s[0], i[0], r[0] = s0, i0, r0
    
    h = t[1] - t[0]  # 今回は 1 週刻み
    for k in range(len(t)-1):
        ds, di, dr = sir_rhs_frac(s[k], i[k], r[k], beta, gamma)
        s[k+1] = s[k] + h*ds
        i[k+1] = i[k] + h*di
        r[k+1] = r[k] + h*dr

        # 数値誤差の安全策
        s[k+1] = max(s[k+1], 0.0)
        i[k+1] = max(i[k+1], 0.0)
        r[k+1] = max(r[k+1], 0.0)

        # 3つが合計1からズレるので軽く正規化（任意）
        total = s[k+1] + i[k+1] + r[k+1]
        if total > 0:
            s[k+1] /= total
            i[k+1] /= total
            r[k+1] /= total

    return s, i, r
```

最適化によってパラメタの値を推定する．
```python
from scipy.optimize import least_squares

# 週インデックス
t = np.arange(len(y) + 1, dtype=float)  
# ↑ 週次新規を作るために s_k - s_{k+1} を使うので、1点多く解を持つ

def predict_weekly_cases(params):
    beta, gamma, i0, C = params
    i0 = float(i0)
    s0 = 1.0 - i0
    r0 = 0.0
    
    s, i, r = euler_sir_frac(s0, i0, r0, beta, gamma, t)
    
    # 週次新規（fraction） = s_k - s_{k+1}
    new_frac = s[:-1] - s[1:]
    y_hat = C * new_frac
    return y_hat, s, i, r

def residuals(params):
    y_hat, _, _, _ = predict_weekly_cases(params)
    return y_hat - y

# 初期値（雑でOK。後で改善できる）
beta0  = 0.8
gamma0 = 0.3
i00    = 1e-4
C0     = y.max() / 0.02  # 適当なスケール初期値

x0 = np.array([beta0, gamma0, i00, C0])

# 境界（不合理な値を防ぐ）
lb = np.array([1e-6, 1e-6, 1e-10, 1.0])
ub = np.array([10.0, 10.0, 0.2, 1e9])

res = least_squares(residuals, x0, bounds=(lb, ub), max_nfev=20000)

beta_hat, gamma_hat, i0_hat, C_hat = res.x
print("=== estimated parameters ===")
print(f"beta  = {beta_hat:.6f}")
print(f"gamma = {gamma_hat:.6f}")
print(f"i0    = {i0_hat:.6e}")
print(f"C     = {C_hat:.6f}")
print(f"R0 = beta/gamma = {beta_hat/gamma_hat:.3f}")
```

