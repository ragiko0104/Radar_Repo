# 微动识别算法 3–9 步的公式与 LaTeX 源

对应 `render_isar_mode_identification_flowchart` 的第 3–9 步，全部与 `functions/` 下的实现一致。

**MathType 导入方式**：MathType 7 → Preferences → Workspace Preferences → 勾选 "Allow TeX language entry from the keyboard"；然后在 MathType 窗口里粘贴下面的 LaTeX 源码（每段一条公式），回 Word 即为可编辑公式。

符号约定：`N` 为慢时间采样点数（锥尖 N = 50，四点联合 N = 200）；`n` 同义；`T_d` 为驻留时长（32 ms）；`f_PR` = 1 kHz；`λ` 为波长（14.8 GHz 对应 0.02027 m）；`t̄` 为 t 的均值。

---

## 3 前向模型与观测算子

### 3.1 姿态链

$$
C(t)=C_0\,R_z(\beta)\,R_y(\theta_p)\,R_z(\alpha)\,R_y(\theta_s)
$$

$$
\beta(t)=-A\cos\left(2\pi f_w t+\psi_w\right),\quad
\theta_p(t)=2\pi f_p t+\psi_p,\quad
\theta_s(t)=2\pi f_s t+\psi_s
$$

```latex
C(t)=C_{0}R_{z}(\beta)R_{y}(\theta_{p})R_{z}(\alpha)R_{y}(\theta_{s})

\beta(t)=-A\cos(2\pi f_{w}t+\psi_{w}),\quad \theta_{p}(t)=2\pi f_{p}t+\psi_{p},\quad \theta_{s}(t)=2\pi f_{s}t+\psi_{s}
```

### 3.2 视线投影

本体坐标 $\mathbf{b}$ 的点经自转、半锥角、进动后投影到视线方向：

$$
\mathbf{u}=R_z(\alpha)R_y(\theta_s)\mathbf{b},\qquad
\mathbf{v}=R_y(\theta_p)\mathbf{u}
$$

$$
R(t)=\frac{1}{\sqrt{2}}\Big[v_1\big(\sin\beta-\cos\beta\big)+v_2\big(\sin\beta+\cos\beta\big)\Big]
$$

```latex
\mathbf{u}=R_{z}(\alpha)R_{y}(\theta_{s})\mathbf{b},\qquad \mathbf{v}=R_{y}(\theta_{p})\mathbf{u}

R(t)=\frac{1}{\sqrt{2}}\left[v_{1}\left(\sin\beta-\cos\beta\right)+v_{2}\left(\sin\beta+\cos\beta\right)\right]
```

锥尖 $\mathbf{b}=(0,-L/2,0)$ 时退化为（$f_s$ 消失）：

$$
R_{\mathrm{tip}}(t)=\frac{L/2}{\sqrt{2}}\Big[\sin\alpha\cos\theta_p\big(\sin\beta-\cos\beta\big)-\cos\alpha\big(\sin\beta+\cos\beta\big)\Big]
$$

```latex
R_{\mathrm{tip}}(t)=\frac{L/2}{\sqrt{2}}\left[\sin\alpha\cos\theta_{p}\left(\sin\beta-\cos\beta\right)-\cos\alpha\left(\sin\beta+\cos\beta\right)\right]
```

### 3.3 观测算子：驻留平均 + 多普勒折叠

$$
f(t)=\frac{2}{\lambda}\cdot\frac{1}{T_d}\int_{-T_d/2}^{T_d/2}\dot R(t+\tau)\,d\tau
$$

$$
\hat f(t)=\frac{2}{\lambda}\cdot\frac{1}{M}\sum_{m=1}^{M}\dot R\!\left(t+\tau_m\right),\qquad M=8
$$

$$
W(x)=x-f_{\mathrm{PR}}\left\lfloor x/f_{\mathrm{PR}}+1/2\right\rfloor
$$

```latex
f(t)=\frac{2}{\lambda}\cdot\frac{1}{T_{d}}\int_{-T_{d}/2}^{T_{d}/2}\dot{R}(t+\tau)\,d\tau

\hat{f}(t)=\frac{2}{\lambda}\cdot\frac{1}{M}\sum_{m=1}^{M}\dot{R}\left(t+\tau_{m}\right),\quad M=8

W(x)=x-f_{\mathrm{PR}}\left\lfloor x/f_{\mathrm{PR}}+1/2\right\rfloor
```

### 3.4 冗余项

$$
\bar f(t;\theta)=W\!\left(\hat f(t)+c_0+c_1(t-\bar t)\right)
$$

```latex
\bar f(t;\theta)=W\!\left(\hat f(t)+c_0+c_1(t-\bar t)\right)
```

---

## 4 假设库与代价函数

$$
\theta=\left[f_p,\ \alpha,\ A,\ f_w,\ \psi_p,\ \psi_w,\ c_0,\ c_1\right]^{\mathsf T}
$$

```
\theta=\left[f_{p},\ \alpha,\ A,\ f_{w},\ \psi_{p},\ \psi_{w},\ c_{0},\ c_{1}\right]^{\mathsf T}
```

$$
\mathrm{RSS}(\theta)=\sum_{m=1}^{N}W\!\left(y_m-\bar f(t_m;\theta)\right)^2
$$

```latex
\mathrm{RSS}(\theta)=\sum_{m=1}^{N}W\!\left(y_m-\bar f(t_m;\theta)\right)^2
```

四个假设（未列出的运动参数被钉为零）：

$$
\mathcal{H}_0:\{c_0,c_1\},\quad
\mathcal{H}_1:\{f_p,\alpha,\psi_p,c_0,c_1\},\quad
\mathcal{H}_2:\{A,f_w,\psi_w,c_0,c_1\},\quad
\mathcal{H}_3:\theta
$$

$$
k(\mathcal H_0)=2,\quad k(\mathcal H_1)=k(\mathcal H_2)=5,\quad k(\mathcal H_3)=8
$$

```latex
\mathcal{H}_{0}:\{c_{0},c_{1}\},\quad \mathcal{H}_{1}:\{f_{p},\alpha,\psi_{p},c_{0},c_{1}\},\quad \mathcal{H}_{2}:\{A,f_{w},\psi_{w},c_{0},c_{1}\},\quad \mathcal{H}_{3}:\theta

k(\mathcal{H}_{0})=2,\quad k(\mathcal{H}_{1})=k(\mathcal{H}_{2})=5,\quad k(\mathcal{H}_{3})=8
```

信息准则（仅作诊断，不参与判决）：

$$
\mathrm{BIC}=N\ln\frac{\mathrm{RSS}}{N}+k\ln N
$$

```latex
\mathrm{BIC}=N\ln\frac{\mathrm{RSS}}{N}+k\ln N
```

---

## 5 选种与双排序并集

网格种子集：

$$
\Theta_{\mathrm{grid}}=\left\{0.4i\right\}_{i=1}^{12}\times\left\{0.4+0.2j\right\}_{j=0}^{22}\times\{2^{\circ},6^{\circ}\}\times\{3^{\circ}\}
$$

$$
\left|\Theta_{\mathrm{grid}}\right|=12\times23\times2=552
$$

```latex
\Theta_{\mathrm{grid}}=\{0.4i\}_{i=1}^{12}\times\{0.4+0.2j\}_{j=0}^{22}\times\{2^{\circ},6^{\circ}\}\times\{3^{\circ}\}

\left|\Theta_{\mathrm{grid}}\right|=12\times23\times2=552
```

谱峰种子：取慢时间谱幅值最大的 3 个 FFT 频点（可能相邻）作 $f_p$，与 $f_w\in\{0.5,1,1.5,2,3,4.5\}$、$\alpha_0\in\{1^{\circ},3^{\circ},5^{\circ},8^{\circ}\}$ 组合，得 72 个，合计 624。

盆地键（同一盆地去重）：

$$
\kappa(\theta)=\frac{1}{20}\,\mathrm{round}\Big(20\left[f_p,\ f_w,\ \alpha\right]^{\mathsf T}\Big)
$$

```latex
\kappa(\theta)=\frac{1}{20}\mathrm{round}\left(20\left[f_{p},\ f_{w},\ \alpha\right]^{\mathsf T}\right)
```

候选盆地集：两套排序各取前 12 个不同键

$$
\mathcal{P}=\mathcal{P}_{\text{short}}\cup\mathcal{P}_{\text{raw}},\qquad
\mathcal{P}_{\text{short}}=\operatorname*{arg\,top12}_{\theta\in\Theta_0}\ \mathrm{RSS}_{\text{50 NM}}(\theta),\quad
\mathcal{P}_{\text{raw}}=\operatorname*{arg\,top12}_{\theta\in\Theta_0}\ \mathrm{RSS}(\theta)
$$

$$
\left|\mathcal P\right|\le 24
$$

```latex
\mathcal{P}=\mathcal{P}_{\mathrm{short}}\cup\mathcal{P}_{\mathrm{raw}},\qquad \left|\mathcal{P}\right|\leq 24
```

---

## 6 有界 Nelder–Mead 精修

H1 在每个候选频率上拟合正弦、余弦、截距和斜率，按 8 点驻留平均换算角度，并按实际观测残差保留候选。H1 的相位在精修中可跨越 ±π，输出时归一化到主值。H3 使用下述有界单纯形。

归一化坐标与尺度：

$$
u_i=\frac{\theta_i}{s_i},\qquad
\mathbf{s}=\left[1,\ 5,\ 5,\ 1,\ 1,\ 1,\ 5,\ 1\right]
$$

```latex
u_{i}=\frac{\theta_{i}}{s_{i}},\qquad \mathbf{s}=\left[1,\ 5,\ 5,\ 1,\ 1,\ 1,\ 5,\ 1\right]
```

初始单纯形（8 个活动参数时为 9 个顶点）：

$$
\mathbf{x}^{(i)}=\mathbf{x}^{(0)}+\delta_i\,\mathbf{e}_i,\qquad
\boldsymbol{\delta}=\left[0.3,\ 0.3,\ 0.3,\ 0.3,\ 1.2,\ 1.2,\ 1.0,\ 0.3\right]
$$

```latex
\mathbf{x}^{(i)}=\mathbf{x}^{(0)}+\delta_{i}\mathbf{e}_{i},\qquad \boldsymbol{\delta}=\left[0.3,\ 0.3,\ 0.3,\ 0.3,\ 1.2,\ 1.2,\ 1.0,\ 0.3\right]
```

四个标准系数（反射 / 扩张 / 收缩 / 缩边）：

$$
\rho=1,\quad \chi=2,\quad \gamma=\tfrac{1}{2},\quad \sigma=\tfrac{1}{2}
$$

```latex
\rho=1,\quad \chi=2,\quad \gamma=\frac{1}{2},\quad \sigma=\frac{1}{2}
```

反射点与扩张点：

$$
\mathbf{x}_r=\bar{\mathbf{x}}+\rho\left(\bar{\mathbf{x}}-\mathbf{x}^{(n+1)}\right),\qquad
\bar{\mathbf{x}}=\frac{1}{n}\sum_{i=1}^{n}\mathbf{x}^{(i)}
$$

$$
\mathbf{x}_e=\bar{\mathbf{x}}+\chi\left(\bar{\mathbf{x}}-\mathbf{x}^{(n+1)}\right)
$$

```latex
\mathbf{x}_{r}=\bar{\mathbf{x}}+\rho\left(\bar{\mathbf{x}}-\mathbf{x}^{(n+1)}\right),\qquad \bar{\mathbf{x}}=\frac{1}{n}\sum_{i=1}^{n}\mathbf{x}^{(i)}

\mathbf{x}_{e}=\bar{\mathbf{x}}+\chi\left(\bar{\mathbf{x}}-\mathbf{x}^{(n+1)}\right)
```

每次移动后做盒约束投影：

$$
\theta_i\leftarrow\min\left(\max\left(\theta_i,\ l_i\right),\ u_i\right)
$$

```latex
\theta_{i}\leftarrow\min\left(\max\left(\theta_{i},\ l_{i}\right),\ u_{i}\right)
```

预算：每个盆地 2500 次评估，取最优后再从该点跑 3000 次；非活动参数钉零（$\theta_i=0,\ i\in\mathcal I_{\text{null}}$）。

---

## 7 F 阶梯判决

四个比较各增加三个参数。N≤8 时入口返回 `RadarProject:TipMode:InsufficientSamples`。

$$
F_1=\frac{\left(\mathrm{RSS}_{\mathcal H_0}-\mathrm{RSS}_{\mathcal H_1}\right)/3}{\mathrm{RSS}_{\mathcal H_1}/(N-5)}
$$

$$
F_2=\frac{\left(\mathrm{RSS}_{\mathcal H_0}-\mathrm{RSS}_{\mathcal H_2}\right)/3}{\mathrm{RSS}_{\mathcal H_2}/(N-5)}
$$

$$
F_3=\frac{\left(\mathrm{RSS}_{\mathcal H_1}-\mathrm{RSS}_{\mathcal H_3}\right)/3}{\mathrm{RSS}_{\mathcal H_3}/(N-8)}
$$

$$
F_4=\frac{(\mathrm{RSS}_{\mathcal H_2}-\mathrm{RSS}_{\mathcal H_3})/3}{\mathrm{RSS}_{\mathcal H_3}/(N-8)}
$$

```latex
F_{1}=\frac{\left(\mathrm{RSS}_{\mathcal{H}_{0}}-\mathrm{RSS}_{\mathcal{H}_{1}}\right)/3}{\mathrm{RSS}_{\mathcal{H}_{1}}/(N-5)}

F_{2}=\frac{\left(\mathrm{RSS}_{\mathcal{H}_{0}}-\mathrm{RSS}_{\mathcal{H}_{2}}\right)/3}{\mathrm{RSS}_{\mathcal{H}_{2}}/(N-5)}

F_{3}=\frac{\left(\mathrm{RSS}_{\mathcal{H}_{1}}-\mathrm{RSS}_{\mathcal{H}_{3}}\right)/3}{\mathrm{RSS}_{\mathcal{H}_{3}}/(N-8)}

F_4=\frac{(\mathrm{RSS}_{\mathcal H_2}-\mathrm{RSS}_{\mathcal H_3})/3}{\mathrm{RSS}_{\mathcal H_3}/(N-8)}
```

等价写法（$F$ 与残差比的关系，$N=50$ 时 $F_3=14\left(\mathrm{RSS}_{\mathcal H_1}/\mathrm{RSS}_{\mathcal H_3}-1\right)$）：

$$
F=\frac{N-k_{\text{big}}}{\Delta k}\left(\frac{\mathrm{RSS}_{\text{small}}}{\mathrm{RSS}_{\text{big}}}-1\right)
$$

```latex
F=\frac{N-k_{\mathrm{big}}}{\Delta k}\left(\frac{\mathrm{RSS}_{\mathrm{small}}}{\mathrm{RSS}_{\mathrm{big}}}-1\right)
```

阈值 τ 为经验值。零幅值时频率和相位不可辨识，未标定这些统计量的零假设分布。

$$
G=(N\ge18)\land(F_3>\tau)\land(F_4\le\tau),\qquad \tau=10
$$

判决按顺序执行：

1. 若 G 成立，选 H2。
2. 否则若 F1 > τ，F3 超过 H1 升级门限时选 H3，其余选 H1。升级门限默认等于 τ，可用 `SwingPromotionF` 单独设置。
3. 否则若 F2 > τ，选 H2；其余选 H0。

```latex
G=(N\ge18)\land(F_3>\tau)\land(F_4\le\tau),\qquad \tau=10
```

---

## 8 自旋频率估计

根部点的多普勒量级（$r_b$ 为根部半径 0.3 m）：

$$
C=\frac{2}{\lambda}\cdot\frac{r_b}{\sqrt 2}\cdot 2\pi=131.5\ \mathrm{Hz/Hz}
$$

$$
\left|f_{\text{base}}\right|\approx C\left(f_s+f_p\right)
$$

```latex
C=\frac{2}{\lambda}\cdot\frac{r_{b}}{\sqrt{2}}\cdot 2\pi=131.5\ \mathrm{Hz/Hz}

\left|f_{\mathrm{base}}\right|\approx C\left(f_{s}+f_{p}\right)
```

由主分量近似得到的搜索上界（代码将上界截在 0–5 Hz）。实际折叠量还取决于角度和姿态。

$$
C\left(f_s+f_p\right)\le\frac{f_{\mathrm{PR}}}{2}
\quad\Longrightarrow\quad
f_s\le\frac{f_{\mathrm{PR}}}{2C}-f_p=3.80-f_p
$$

```latex
C\left(f_{s}+f_{p}\right)\leq\frac{f_{\mathrm{PR}}}{2}\quad\Longrightarrow\quad f_{s}\leq\frac{f_{\mathrm{PR}}}{2C}-f_{p}=3.80-f_{p}
```

联合目标（固定第一级结果，逐点仿射冗余被消去）：

$$
\mathrm{RSS}_s\left(f_s,\psi_s,\alpha,\psi_p\right)=\sum_{k=1}^{4}\left\|\mathbf{d}_k-A_2\left(A_2^{\dagger}\mathbf{d}_k\right)\right\|^2
$$

$$
\mathbf{d}_k=W(\mathbf{y}_k-\mathbf{M}_k),\qquad
A_2=\left[\mathbf{1},\ \mathbf{t}-\bar t\,\mathbf{1}\right]
$$

```latex
\mathrm{RSS}_{s}\left(f_{s},\psi_{s},\alpha,\psi_{p}\right)=\sum_{k=1}^{4}\left\|\mathbf{d}_{k}-A_{2}\left(A_{2}^{\dagger}\mathbf{d}_{k}\right)\right\|^{2}

\mathbf{d}_k=W(\mathbf{y}_k-\mathbf{M}_k),\qquad A_2=[\mathbf{1},\mathbf{t}-\bar t\,\mathbf{1}]
```

无自旋模型固定 fs=0，另外三个非线性参数可重拟合。两模型各自逐轨估计截距和斜率。四轨满秩时，大模型的名义残差自由度为 N−4−2×4=N−12，参数数目差为 1。

$$
F_s=\frac{\mathrm{RSS}_{s,0}-\mathrm{RSS}_{s,1}}{\mathrm{RSS}_{s,1}/(N-12)}
$$

```latex
F_{s}=\frac{\mathrm{RSS}_{s,0}-\mathrm{RSS}_{s,1}}{\mathrm{RSS}_{s,1}/(N-12)}
```

识别条件（三者同时成立）：

$$
F_s>10\ \wedge\ \hat f_s>0.2\ \mathrm{Hz}\ \wedge\ \sqrt{\mathrm{RSS}_{s,1}/N}<10\ \mathrm{Hz}
$$

```latex
F_{s}>10\ \wedge\ \hat{f}_{s}>0.2\ \mathrm{Hz}\ \wedge\ \sqrt{\mathrm{RSS}_{s,1}/N}<10\ \mathrm{Hz}
```

---

## 9 输出与不确定度

残差方差与协方差（$J$ 为残差对活动参数的有限差分雅可比）：

$$
\hat\sigma^2=\frac{\mathrm{RSS}}{N-k},\qquad
\mathrm{Cov}\left(\hat\theta\right)=\hat\sigma^2\left(J^{\mathsf T}J\right)^{\dagger}
$$

$$
\mathrm{sd}\left(\hat\theta_i\right)=\sqrt{\left[\mathrm{Cov}\right]_{ii}},\qquad
J_{mi}=\frac{\partial r_m}{\partial\theta_i}\approx\frac{r_m\left(\theta_i+\Delta_i\right)-r_m\left(\theta_i\right)}{\Delta_i}
$$

```latex
\hat{\sigma}^{2}=\frac{\mathrm{RSS}}{N-k},\qquad \mathrm{Cov}\left(\hat{\theta}\right)=\hat{\sigma}^{2}\left(J^{\mathsf T}J\right)^{\dagger}

\mathrm{sd}\left(\hat{\theta}_{i}\right)=\sqrt{\left[\mathrm{Cov}\right]_{ii}},\qquad J_{mi}=\frac{\partial r_{m}}{\partial\theta_{i}}\approx\frac{r_{m}\left(\theta_{i}+\Delta_{i}\right)-r_{m}\left(\theta_{i}\right)}{\Delta_{i}}
```

锥尖返回最终模型、参数、局部标准误和 F 阶梯。自旋表返回 `estSpinHz`、`spinDetected` 和 `jointRmsHz`。`aliased` 来自场景真值，仅用于评估；通过检出门限的样本也可能折叠。`fw≈fp` 时仍返回普通模式和参数，使用者须结合退化情形解释。

---

## 附：第 3–4 步的补充——旋转矩阵、解析导数与变量表

### A.1 旋转矩阵定义

$$
R_z(\phi)=\begin{bmatrix}\cos\phi&-\sin\phi&0\\ \sin\phi&\cos\phi&0\\ 0&0&1\end{bmatrix},\qquad
R_y(\phi)=\begin{bmatrix}\cos\phi&0&\sin\phi\\ 0&1&0\\ -\sin\phi&0&\cos\phi\end{bmatrix}
$$

```latex
R_{z}(\phi)=\begin{bmatrix}\cos\phi&-\sin\phi&0\\ \sin\phi&\cos\phi&0\\ 0&0&1\end{bmatrix},\qquad R_{y}(\phi)=\begin{bmatrix}\cos\phi&0&\sin\phi\\ 0&1&0\\ -\sin\phi&0&\cos\phi\end{bmatrix}
```

### A.2 视线投影的解析导数

$$
\dot{\mathbf{u}}=\dot\theta_s\,R_z(\alpha)\,R_y'(\theta_s)\,\mathbf{b},\qquad
\dot{\mathbf{v}}=\dot\theta_p\,R_y'(\theta_p)\,\mathbf{u}+R_y(\theta_p)\,\dot{\mathbf{u}}
$$

$$
\dot R(t)=\frac{1}{\sqrt2}\Big[\left(\dot v_1-\dot\beta v_2\right)\left(\sin\beta-\cos\beta\right)+\left(\dot v_2+\dot\beta v_1\right)\left(\sin\beta+\cos\beta\right)\Big]
$$

$$
\dot\beta=2\pi f_w A\sin\left(2\pi f_w t+\psi_w\right),\qquad \dot\theta_p=2\pi f_p,\qquad \dot\theta_s=2\pi f_s
$$

```latex
\dot{\mathbf{u}}=\dot{\theta}_{s}R_{z}(\alpha)R_{y}'(\theta_{s})\mathbf{b},\qquad \dot{\mathbf{v}}=\dot{\theta}_{p}R_{y}'(\theta_{p})\mathbf{u}+R_{y}(\theta_{p})\dot{\mathbf{u}}

\dot{R}(t)=\frac{1}{\sqrt{2}}\left[\left(\dot{v}_{1}-\dot{\beta}v_{2}\right)\left(\sin\beta-\cos\beta\right)+\left(\dot{v}_{2}+\dot{\beta}v_{1}\right)\left(\sin\beta+\cos\beta\right)\right]

\dot{\beta}=2\pi f_{w}A\sin\left(2\pi f_{w}t+\psi_{w}\right),\qquad \dot{\theta}_{p}=2\pi f_{p},\qquad \dot{\theta}_{s}=2\pi f_{s}
```

### A.3 变量表

| 符号 | 含义 | 单位 | 取值/来源 |
|---|---|---|---|
| $t$ | 慢时间 | s | 0–5 |
| $t_m$ | 第 $m$ 帧 CPI 中心时刻 | s | 合成扫描采用 $0.0156+0.1(m-1)$；提取轨迹读取实际接收中心时刻 |
| $N$ | 帧数（慢时间采样点数） | — | 50 |
| $\lambda$ | 波长 | m | $c/14.8\,\mathrm{GHz}=0.02027$ |
| $T_d$ | 相参驻留时长 | s | $32\times10^{-3}$ |
| $M$ | 驻留平均的子样本数 | — | 8 |
| $f_{\mathrm{PR}}$ | 脉冲重复频率 | Hz | 1000 |
| $C(t)$ | 本体到视线的姿态矩阵 | — | 由下式给出 |
| $C_0$ | 初始姿态 | — | $R_z(-135^\circ)$ |
| $\beta(t)$ | 摆动角 | rad | $-A\cos(2\pi f_w t+\psi_w)$ |
| $A$ | 摆幅 | rad（界面用度） | 待估，0–25° |
| $f_w$ | 摆动频率 | Hz | 待估，0.3–4.9 |
| $\psi_w$ | 摆动初相 | rad | 待估 |
| $\theta_p(t)$ | 进动角 | rad | $2\pi f_p t+\psi_p$ |
| $f_p$ | 进动频率 | Hz | 待估，0.3–4.9 |
| $\psi_p$ | 进动初相 | rad | 待估 |
| $\theta_s(t)$ | 自转角 | rad | $2\pi f_s t+\psi_s$ |
| $f_s$ | 自转频率 | Hz | 第 8 步单独估 |
| $\psi_s$ | 自转初相 | rad | 第 8 步单独估 |
| $\alpha$ | 进动半锥角 | rad（界面用度） | 待估，0–25° |
| $\mathbf{b}$ | 散射点本体坐标 | m | 锥尖 $(0,-0.75,0)$ |
| $\mathbf{u},\mathbf{v}$ | 中间坐标 | m | 见 A.2 |
| $R(t)$ | 视线投影距离 | m | — |
| $\dot R(t)$ | 距离变化率 | m/s | 解析导数 |
| $f(t)$ | 瞬时多普勒 | Hz | $2\dot R/\lambda$ |
| $f_{\mathrm{PR}}$ 折叠后 $\tilde f$ | 无模糊区间内的多普勒 | Hz | $\pm500$ |
| $c_0$ | 零阶冗余（距离门带来的多普勒偏置） | Hz | 待估，约 5 |
| $c_1$ | 一阶冗余（距离游走耦合） | Hz/s | 待估，约 0.14 |
| $\bar f(t;\theta)$ | 完整模型输出 | Hz | — |
| $y_m$ | 第 $m$ 帧提取的多普勒 | Hz | 数据，第 2 步产出 |
| $r_m$ | 残差 | Hz | $W(y_m-\bar f(t_m;\theta))$ |
| $\theta$ | 待估参数向量 | — | 8 维 |
| $k$ | 自由参数个数 | — | 2/5/5/8 |
