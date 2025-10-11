# 物理学の事例

到達目標
- 指数的減衰を示す現象を数理的に理解する．

準備
1. anacondaを使用し，<span style="color:red">jupyter lab</span>を起動する．
2. `Documents（書類）/mathematical_studies`フォルダをダブルクリックで開き`+`をクリックして新しいファイルを作成する．
3. ファイル名を`2_{学籍番号}_{氏名}.ipynb`に変更する．例：`2_SI25999_香川渓一郎.ipynb`

## 前回の内容

数理モデルは次の3要素から構成される．

- **変数**：対象となる現象の状態・性質・量などを数字やラベルで表したもの．
- **数理構造**：変数が従うルールを数学的に表現したもの．数理モデルの骨格．
- **パラメタ**：数理モデルの表現可能性を担う定数．これを調節することでデータを説明する．

人口増加の数理モデルは次の微分方程式によって表される．

$$
\frac{dN}{dt} = r N(t)
$$

ここに，$N(t)$ は各時刻 $t$ における人口を表す変数であり，$r$ は増加率を表すパラメタで $r > 0$ だった．

## 放射性崩壊

### 歴史

- 1858年：プリュッカー (Julius Plücker) による陰極線の発見
  - 真空管（クルックス管）に電気を通すと，陽極の背後のガラスが緑色に発光する現象が発見された．
  - 陰極から何かしらの放射線が出ていることを示唆しており，これを陰極線と呼んだ．
- 1895年：レントゲン (Wilhelm Conrad Röntgen) による<span style="color:red">X線</span>の発見
  - 真空管にアルミ箔の窓をつけると，陰極線が真空管の外に出ることが知られていた．
  - 真空管の外に出た陰極線は2メートル離れた蛍光板を光らせることを発見した．
  - 真空管を黒い紙で覆い，蛍光板との間に本，木片，ガラスなどを挟んでも蛍光板は光ることから，これは未知の放射線であるとしX線と名付けた．
  - この発見によりレントゲンは1901年の第1回ノーベル物理学賞を受賞する．

![アルミ箔付きクルックス管](/contents/figs/2/Crookes_tube.jpg)

[Author: D-Kuru/Wikimedia Commons, Licence: CC-BY-SA-3.0-AT](https://commons.wikimedia.org/wiki/File:Crookes_tube-in_use-lateral_view-standing_cross_prPNr%C2%B011.jpg)より

- 1896年：ベクレル (Antoine Henri Becquerel) による<span style="color:red">放射能</span>の発見
  - 蛍光物質として知られていた**ウラン化合物**に日光を当てればX線が放射されると考え，それを実験により確認した．
  - しかしながら，日光を当てなくてもX線より強い放射線が放出されていることを偶然発見した．
  - このことから，ウラン化合物が**放射能**（放射線を出す能力）を持っていることを発見する．
  - ベクレルの名前は放射能の強さを表す単位に残る．

### 放射性崩壊の原理

![放射性崩壊のイメージ](/contents/figs/2/01-02-05.png)

[環境省「放射線による健康影響等に関する統一的な基礎資料」](https://www.env.go.jp/chemi/rhm/current/kisoshiryohtml.html)より

![放射性崩壊の種類](/contents/figs/2/01-02-06.png)

[環境省「放射線による健康影響等に関する統一的な基礎資料」](https://www.env.go.jp/chemi/rhm/current/kisoshiryohtml.html)より

ここでは不安定核種が直接安定核種に変化する場合を扱う．

### 数理モデルの構築

- 単位時間あたりに放射性崩壊する確率（**崩壊定数**）を $\lambda > 0$ とし，時刻 $t$ での放射性物質の数を $N(t)$ とおく．
- 現在時刻を $t=0$ として，現在時刻における放射性物質の数が $N_0$ と分かっているものとする．
- このとき，任意の時刻 $t$ での放射性物質の数 $N(t)$ と，微小時間 $\varepsilon$ 後の放射性物質の数 $N(t+\varepsilon)$ の関係は次のように求められる．

$$
N(t+\varepsilon) = N(t) \times (1- \lambda \varepsilon).
$$

これを式変形すれば

$$
\frac{N(t+\varepsilon)-N(t)}{\varepsilon} = - \lambda N(t)．
$$

となる．
ここで $\varepsilon \to 0$ の極限を取れば，次の微分方程式を得る．

$$
\frac{dN}{dt} = - \lambda N(t).
\tag{1}
$$

これを解くと（この関係式を満たす $N(t)$ を求めると）

$$
N(t) = N_0 e^{- \lambda t}.
$$

### 半減期

放射性崩壊によって粒子数が半分になるまでにかかる時間を<span style="color:red">半減期</span>と呼ぶ．
半減期によって人体に及ぼす影響を評価することができる．

半減期を $T$ とするとき，次の関係を満たす $T$ を求めれば良い．

$$
\frac{N_0}{2} = N_0 \ e^{- \lambda T}.
$$

従って

$$
T = \frac{\log 2}{\lambda}
$$

と求まる．

### グラフのプロット

準備
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'
```

変数の準備
```python
# 初期値と崩壊定数を定める
N_0 = 10000
r = 0.01
```

プロット
```python
# FigureとAxesを生成する
fig, ax = plt.subplots(figsize=(5, 3))

# 時間軸を定める
t = np.linspace(0, 1, 100)

# 横軸を時間，縦軸をN(t)としてプロットする
ax.plot(t, N_0 * np.exp(-r * t))
ax.set_xlabel('時刻'); ax.set_ylabel('放射性物質の数')

# Figureを保存する（相対パスを指定）
fig.savefig("./2_decay.png")
```

```{note}
**演習1**

1. トリチウムの崩壊定数は $\lambda = 1.8 \times 10^{-9} \mathrm{[s^{-1}]}$ であり，初期値の粒子数は $5.6 \times 10^{12}$ 個とする．
このとき，半減期は何秒（何年）になるか．

※ $5.6 \times 10^{12}$ 個は世界保健機関 (WHO) による飲料水1リットルの水質の基準値．

2. このときのグラフを作成せよ．
```

## 冷却現象

### 歴史

- 17世紀末：ニュートン (Isaac Newton) による冷却の法則の発見[3]
  - 熱い液体の温度が時間とともにどのように下がるかを測定した．
  - 液体が冷える速さは、周囲との温度差が大きいほど早く，その速さは液体と周囲の温度との**温度差に比例する**ことを見出した．

![ニュートン](/contents/figs/2/GodfreyKneller-IsaacNewton-1689.jpg)

[Public Domain / Wikimedia Commons](https://commons.wikimedia.org/wiki/File:GodfreyKneller-IsaacNewton-1689.jpg)

- この法則は、日常的な冷却現象だけでなく、**熱伝達・気象・生体熱交換**など多くの分野で応用されている。

### 冷却現象の原理

- 伝熱の3つの形態：
  1. **伝導 (conduction)**：固体内での熱の移動
  2. **対流 (convection)**：空気や水の流れによる熱輸送
  3. **放射 (radiation)**：赤外線などの電磁波による熱放出

- ニュートンの冷却則は，**対流による放熱が支配的な場合**に成立する近似式である．実際，ニュートンは風を当てながら温度の低下を調べた．

### 数理モデルの構築

- 時刻 $t$ における物体の温度を $T(t)$，周囲の温度を $T_{\mathrm{env}}$ とする．
- 実験的事実として**温度変化の速さ**は**温度差**に比例することが知られている．
- 現在時刻を $t=0$ として，現在時刻における物体の温度が $T_0 > T_{\mathrm{env}}$ と分かっているものとする．
- 温度の下がり方を決める定数（冷却定数）を $k > 0$ とする．

このとき

$$
\frac{dT}{dt} = -k \bigl(T(t) - T_{\mathrm{env}}\bigr).
\tag{1}
$$

これを解くと（この関係式を満たす $T(t)$ を求めると）

$$
\log|T(t) - T_{\mathrm{env}}| - \log(T_0 - T_{\mathrm{env}}) = -kt.
$$

指数関数をとると、

$$
T(t) = T_{\mathrm{env}} + (T_0 - T_{\mathrm{env}}) e^{-kt}.
$$

### 特徴

- 時間が経つと物体の温度は周囲の温度に漸近する（$\lim_{t\to 0} T(t) = T_{\mathrm{env}}$）．
- 物体と周囲との温度差は指数関数的に減少する：

$$
T(t) - T_{\mathrm{env}} = (T_0 - T_{\mathrm{env}}) e^{-kt}.
$$

- 放射性崩壊と同様，「変化率が現在値に比例」する指数的減衰モデルである．

### グラフのプロット

準備
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Hiragino Sans'
```

変数の準備
```python
# 初期温度，環境温度，冷却定数を設定
T_env = 25     # [℃]
T_0 = 90       # [℃]
k1 = 0.01       # [1/s]
k2 = 0.02       # [1/s]
```

プロット
```python
# FigureとAxesを生成
fig, ax = plt.subplots(figsize=(5,3))

# 時間軸
t = np.linspace(0, 1000, 100)
# 温度の時間変化を計算
T1 = T_env + (T_0 - T_env) * np.exp(-k1 * t)
T2 = T_env + (T_0 - T_env) * np.exp(-k2 * t)

# プロット
ax.plot(t, T1, label=f"温度変化 (k={k1})")
ax.plot(t, T2, label=f"温度変化 (k={k2})")
ax.set_xlabel("時間 [s]")
ax.set_ylabel("温度 [${}^\circ\mathrm{C}$]")
ax.legend()

fig.savefig("./2_cooling.png")
```

```{note}
**演習2**

1. 初期温度 $T_0 = 90[{}^\circ\mathrm{C}]$，環境温度 $T_{\mathrm{env}} = 25[{}^\circ\mathrm{C}]$ とする．1000秒後に温度が $60[{}^\circ\mathrm{C}]$ になったとき，冷却定数 $k$ の値を求めよ．

2. そのときのパラメタ設定を用いて2000秒後の温度を求めよ．

3. 異なる冷却定数 $k$ の値（例：0.008, 0.005）でグラフを重ね，冷却速度の違いを比較せよ．
```

## 補足：指数的減衰との共通性

| 現象      | 方程式                                        | 解                                                          |
| -----    | -------------------------------------------- | ----------------------------------------------------------- |
| 放射性崩壊 | $\frac{dN}{dt} = -\lambda N(t)$              | $N(t) = N_0 e^{-\lambda t}$                                 |
| 冷却現象  | $\frac{dT}{dt} = -k(T - T_{\mathrm{env}})$   | $T(t) = T_{\mathrm{env}} + (T_0 - T_{\mathrm{env}})e^{-kt}$  |

<u>共通点</u>：変化率が「現在値（または現在の温度差）」に比例 → 指数関数的減衰。

<u>違い</u>：冷却現象では環境との「温度差」が減衰する．

## まとめ

本日作成したipynbファイルに追記する形で次の課題に回答し，WebClassの第2回課題から提出せよ．

```{note}
**復習**

(1) 本日学んだことを3つ箇条書きで述べよ．
```

## 参考文献

1. [内閣官房「放射線研究の幕開け ～レントゲンによるX線の発見～」](https://www.kantei.go.jp/saigai/senmonka_g51.html)
2. [環境省「放射線による健康影響等に関する統一的な基礎資料」](https://www.env.go.jp/chemi/rhm/current/kisoshiryohtml.html)
3. “Scala graduum Caloris (A Scale of the Degrees of Heat)”, Philosophical Transactions, No. 270, pp. 824–829, (April 1701). 
4. [円山重直「ニュートンの冷却法則（その 1）」](http://www.wattandedison.com/Maruyama2015.pdf)
