# Fullpath 六段 LFM 的 SNR/SINR 测量精度门

本文说明 `multiuav_rvpath.m` 中六段 LFM 距离/速度（R/V）解算的一致性门限原理、公式推导、候选关联约束和配置方法。

## 1. 设计目的

六段 LFM 链路从六个频段的脉压峰中选出能量最高的三个频段，分别组成两组双频段 R/V 解：

- 第一组使用峰时延 \(\tau_1,\tau_2\)，得到 \(R_{12},V_{12}\)；
- 第二组使用峰时延 \(\tau_1,\tau_3\)，得到 \(R_{13},V_{13}\)；
- 两组解共享第一个频段的峰。

如果三个峰来自同一个目标，两组结果应当在测量误差范围内一致。门限的任务是判断：

\[
R_{12}\approx R_{13},\qquad V_{12}\approx V_{13}.
\]

门限只能判断两组结果在统计意义上是否相容，不能独立判断三个峰是否来自同一个物理目标。因此，完整链路还需要距离顺序约束、峰的一对一使用约束和物理速度硬门限。

## 2. Legacy 固定门限

legacy 优化链使用固定距离门和速度门：

```matlab
Rres = c/min(abs([B1,B2,B3]))/5;
Vres = (c/2/(Tr*fc))*2/5;
```

当前参数为：

```text
B = [5, 7.5, 10] MHz
Tr = 50 us
fc = 12.5 GHz
```

因此：

\[
R_{\mathrm{legacy}}=12\ \mathrm{m},
\qquad
V_{\mathrm{legacy}}=96\ \mathrm{m/s}.
\]

这两个门限是固定分辨尺度的经验比例，不会随当前脉压峰的 SNR、目标间干扰和峰拟合精度变化。

## 3. 总体处理流程

新的 fullpath R/V 候选判定流程为：

```text
脉压峰功率与背景功率
        ↓
估计每个峰的有效 SNR/SINR
        ↓
计算每个峰的时延标准差 sigmaTau
        ↓
通过 Cal_rv_mat 雅可比传播为 sigmaDeltaR、sigmaDeltaV
        ↓
校正确定性跨频段偏差 biasR、biasV
        ↓
执行 N-sigma 距离/速度一致性门
        ↓
执行跨频段峰距离顺序约束
        ↓
执行三个频段峰的一对一使用约束
        ↓
执行 ±20 Mach 物理速度硬门限
```

## 4. 峰值 SNR/SINR 估计

### 4.1 基本定义

对第 \(i\) 个脉压峰，线性功率域的有效 SNR/SINR 定义为：

\[
\gamma_i
=
\max\left(
\frac{P_{\mathrm{peak},i}}
     {P_{\mathrm{noise},i}}
-1,
\gamma_{\min}
\right).
\]

同时设置上限：

\[
\gamma_i\leq \gamma_{\max},
\]

避免无噪仿真或数值零噪声导致无限 SNR 和零宽门限。

默认配置为：

```matlab
snrFloorDb   = 0;
snrCeilingDb = 60;
```

### 4.2 `localCfar`：默认有效 SINR

默认噪声估计方式为：

```matlab
noiseEstimator = 'localCfar';
```

CFAR 阈值满足：

\[
T_{\mathrm{CFAR}}=\alpha P_{\mathrm{training}},
\]

所以峰附近的训练区背景功率为：

\[
P_{\mathrm{noise},i}
=
\frac{T_{\mathrm{CFAR},i}}{\alpha}.
\]

这个背景中不仅包含接收噪声，还可能包含：

- 其他目标的脉压旁瓣；
- 多散射点相干干扰；
- 邻近目标泄漏；
- 数值链背景。

因此，`localCfar` 得到的更准确含义是有效 SINR。它通常比纯热噪声 SNR 更能反映多目标环境下峰位置拟合的实际精度。

### 4.3 `robustOffTarget`：离目标稳健噪声底

可选配置为：

```matlab
noiseEstimator = 'robustOffTarget';
```

程序排除所有已检测峰及其主瓣/保护区域，再使用剩余脉压功率样本的中位数估计背景。对于复高斯噪声，功率服从指数分布，其均值估计为：

\[
P_{\mathrm{noise}}
\approx
\frac{\operatorname{median}(P)}{\ln 2}.
\]

这种方式更接近纯噪声底，但不会把确定性旁瓣当成干扰。在强多目标场景中，它可能高估峰拟合精度，使门限过窄。

## 5. 从 SINR 推导时延测量精度

### 5.1 信号均方根带宽

设基带信号频谱为 \(S(f)\)，其均方根带宽定义为：

\[
\beta_{\mathrm{rms}}^2
=
\frac{
\int (f-\bar f)^2|S(f)|^2\,\mathrm{d}f
}{
\int |S(f)|^2\,\mathrm{d}f
}.
\]

对于带宽为 \(B\) 的理想矩形 LFM 频谱：

\[
\beta_{\mathrm{rms}}=\frac{B}{\sqrt{12}}.
\]

### 5.2 时延估计的 Cramér–Rao 下界

在复高斯白噪声、无偏时延估计条件下，时延参数的 Fisher 信息为：

\[
\mathcal I_{\tau}
=
8\pi^2\beta_{\mathrm{rms}}^2\gamma.
\]

所以时延估计方差满足：

\[
\operatorname{var}(\hat\tau)
\geq
\frac{1}
{8\pi^2\beta_{\mathrm{rms}}^2\gamma}.
\]

对应的标准差为：

\[
\sigma_{\tau,\mathrm{CRLB}}
=
\frac{1}
{2\sqrt{2}\pi\beta_{\mathrm{rms}}\sqrt{\gamma}}.
\]

因此：

\[
\sigma_\tau
\propto
\frac{1}{B\sqrt{\gamma}}.
\]

这表示：

- 带宽增大，时延精度提高；
- SINR 增大，时延精度提高；
- SINR 提高 6 dB，即线性 SINR 增大约 4 倍，随机时延标准差约减半。

### 5.3 峰拟合与数值模型精度下限

CRLB 不包含有限采样、三点峰拟合、在线频谱近似和模型失配。程序增加一个以采样 bin 为单位的标准差下限：

\[
\sigma_{\tau,\mathrm{floor}}
=
\frac{\sigma_{\mathrm{floor,bin}}}{f_s}.
\]

最终时延标准差为：

\[
\sigma_{\tau}
=
\sqrt{
\sigma_{\tau,\mathrm{CRLB}}^2
+
\sigma_{\tau,\mathrm{floor}}^2
}.
\]

当前默认值为：

```matlab
delayStdFloorBins = 0.05;
```

## 6. 双频段距离/速度解算

对频段 \(a,b\)，定义：

\[
q_a=\frac{k_a}{f_a},
\qquad
q_b=\frac{k_b}{f_b},
\]

其中 \(k\) 为 LFM 调频斜率，\(f\) 为对应频段载频。

`Cal_rv_mat` 的距离解为：

\[
R_{ab}
=
\frac{c}{2}
\frac{q_a\tau_a-q_b\tau_b}
{q_a-q_b}.
\]

速度解为：

\[
V_{ab}
=
\frac{c}{2}
\frac{\tau_a-\tau_b}
{f_a/k_a-f_b/k_b}.
\]

使用三个频段得到两组解：

\[
(R_{12},V_{12}),
\qquad
(R_{13},V_{13}).
\]

## 7. `Cal_rv_mat` 误差的解析传播

### 7.1 距离解偏导

定义：

\[
A_{ab,a}
=
\frac{c}{2}\frac{q_a}{q_a-q_b},
\]

\[
A_{ab,b}
=
-\frac{c}{2}\frac{q_b}{q_a-q_b}.
\]

则：

\[
R_{ab}=A_{ab,a}\tau_a+A_{ab,b}\tau_b.
\]

距离差定义为：

\[
\Delta R=R_{12}-R_{13}.
\]

它对三个时延的雅可比为：

\[
J_R
=
\begin{bmatrix}
A_{12,1}-A_{13,1} & A_{12,2} & -A_{13,3}
\end{bmatrix}.
\]

### 7.2 速度解偏导

定义：

\[
G_{ab}
=
\frac{c/2}{f_a/k_a-f_b/k_b}.
\]

则：

\[
V_{ab}=G_{ab}(\tau_a-\tau_b).
\]

速度差定义为：

\[
\Delta V=V_{12}-V_{13}.
\]

它对三个时延的雅可比为：

\[
J_V
=
\begin{bmatrix}
G_{12}-G_{13} & -G_{12} & G_{13}
\end{bmatrix}.
\]

### 7.3 标准差传播

假设三个频段的峰时延误差近似独立：

\[
\Sigma_\tau
=
\operatorname{diag}
\left(
\sigma_{\tau_1}^2,
\sigma_{\tau_2}^2,
\sigma_{\tau_3}^2
\right).
\]

距离差和速度差的方差为：

\[
\sigma_{\Delta R}^2
=
J_R\Sigma_\tau J_R^{\mathrm T},
\]

\[
\sigma_{\Delta V}^2
=
J_V\Sigma_\tau J_V^{\mathrm T}.
\]

展开后：

\[
\sigma_{\Delta R}
=
\sqrt{
\sum_{i=1}^{3}
(J_{R,i}\sigma_{\tau_i})^2
},
\]

\[
\sigma_{\Delta V}
=
\sqrt{
\sum_{i=1}^{3}
(J_{V,i}\sigma_{\tau_i})^2
}.
\]

因此，每一组三峰候选都有由其三个峰 SINR 决定的独立距离门和速度门。

## 8. 确定性跨频段偏差校正

CRLB 描述随机噪声误差，不包括以下确定性误差：

- 在线 LFM 频谱的频率 bin 取整；
- 不同带宽下三点峰拟合的系统偏差；
- 不同频段的固定时延参考差；
- 多散射点相干叠加造成的可重复峰偏移。

程序默认使用：

```matlab
biasCorrection = 'rankMedian';
```

先将三个频段的峰按距离从近到远排列。对相同序号的峰计算：

\[
dR_i=R_{12}(i,i)-R_{13}(i,i),
\]

\[
dV_i=V_{12}(i,i)-V_{13}(i,i).
\]

再使用稳健中位数估计公共偏差：

\[
b_R=\operatorname{median}(dR_i),
\]

\[
b_V=\operatorname{median}(dV_i).
\]

门限针对校正后的残差：

\[
\varepsilon_R
=
|R_{12}-R_{13}-b_R|,
\]

\[
\varepsilon_V
=
|V_{12}-V_{13}-b_V|.
\]

该校正只消除跨频段公共偏差，不使用距离或速度理论真值。

## 9. N-sigma 一致性门

配置：

```matlab
sigmaMultiplier = 2;
```

候选通过条件为：

\[
\varepsilon_R
<
N_\sigma\sigma_{\Delta R},
\]

\[
\varepsilon_V
<
N_\sigma\sigma_{\Delta V}.
\]

当前 \(N_\sigma=2\)。若单变量误差近似高斯，双侧 \(2\sigma\) 区间覆盖概率约为 95.45%。距离误差和速度误差来自相同的三个时延，并非严格独立，所以不能简单把两个 95.45% 相乘作为精确联合通过概率。

候选归一化代价为：

\[
D
=
\sqrt{
\left(
\frac{\varepsilon_R}
{N_\sigma\sigma_{\Delta R}}
\right)^2
+
\left(
\frac{\varepsilon_V}
{N_\sigma\sigma_{\Delta V}}
\right)^2
}.
\]

同一个公共峰存在多个相容候选时，优先选择 \(D\) 最小的第三频段峰。

## 10. 为什么不能只缩小门限

当前十目标四点圆周仿真表明：

- legacy 固定门产生 16 个候选，其中 6 个为组合假目标；
- 把门收紧到纯噪声 CRLB 对应的约 `0.713 m / 11.33 m/s` 后，只剩 3 个候选；
- 这 3 个候选仍然是组合假目标；
- 一些组合假目标的两组解算残差比真实目标更小。

因此，统计一致性并不等同于物理同源性。只缩门限可能优先拒绝受跨频段模型偏差影响的真实目标，而保留偶然高度一致的组合假目标。

## 11. 距离顺序和一对一约束

### 11.1 距离顺序约束

当前默认：

```matlab
maxPeakRankDifference = 0;
```

表示三个频段只允许相同距离排序序号的峰参与组合：

\[
i_1=i_2=i_3.
\]

当实测数据可能发生漏峰时，可以设置：

```matlab
maxPeakRankDifference = 1;
```

允许相邻序号峰匹配。设置为 `Inf` 可以关闭距离排序约束。

### 11.2 峰的一对一使用

当前默认：

```matlab
enforceOneToOnePeakUse = true;
```

所有通过一致性门的候选按照归一化代价 \(D\) 从小到大排序。候选被接受后，它使用的三个频段峰都会被标记为已占用，不能再次参与另一个目标组合。

该约束防止一个强峰与多个其他频段峰重复组合成多个假目标。

## 12. 物理速度硬门限

统计门和关联约束之后，仍保留物理速度门：

\[
|V|\leq 20\,\mathrm{Mach}.
\]

以 `Mach=340 m/s` 计：

\[
|V|\leq6800\ \mathrm{m/s}.
\]

该门用于剔除物理上不可能的极端速度解，不用于替代跨频段一致性判断。

## 13. 当前默认配置

配置位于 `code/multiuav_iasr_fullpath.m`：

```matlab
cfg_RVpath.rvConsistency = struct( ...
    'mode','snrPrecision', ...
    'sigmaMultiplier',2, ...
    'snrFloorDb',0, ...
    'snrCeilingDb',60, ...
    'delayStdFloorBins',0.05, ...
    'noiseEstimator','localCfar', ...
    'biasCorrection','rankMedian', ...
    'maxPeakRankDifference',0, ...
    'enforceOneToOnePeakUse',true);
```

兼容 legacy 固定门可以设置：

```matlab
cfg_RVpath.rvConsistency.mode = 'legacyResolution';
```

若需要尽量复现旧的全组合行为，还应设置：

```matlab
cfg_RVpath.rvConsistency.biasCorrection = 'none';
cfg_RVpath.rvConsistency.maxPeakRankDifference = Inf;
cfg_RVpath.rvConsistency.enforceOneToOnePeakUse = false;
```

## 14. 当前仿真结果

使用十目标、四点圆周散射模型、无附加噪声的当前 fullpath 参数：

```text
有效 SINR              14.72 至 19.62 dB
候选距离门             8.07 至 9.83 m
候选速度门             128.38 至 156.32 m/s
最终候选数             10
额外候选数             0
漏检数                 0
最大原始距离相对误差   0.2533%
最大原始速度相对误差   0.7188%
```

速度门大于 legacy 的 `96 m/s` 并不意味着门限失效。速度由多个频段的时延差反演；当对应雅可比增益较大时，很小的时延误差会被放大为较大的速度标准差。当前组合假目标主要由距离顺序和一对一约束消除，而不是通过盲目缩窄速度门消除。

默认 `quad_small` UAV 场景的当前结果为：

```text
最终候选数             10
额外候选数             0
漏检数                 0
最大原始距离相对误差   0.504%
最大原始速度相对误差   2.63%
```

## 15. 输出诊断字段

主要诊断量保存在 `rvMeasurement`：

```matlab
rvMeasurement.parameters
rvMeasurement.selectedPeakSnrDb
rvMeasurement.selectedPeakNoisePower
rvMeasurement.selectedPeakLocalCfarNoisePower
rvMeasurement.selectedPeakRobustOffTargetNoisePower
rvMeasurement.selectedPeakDelayStdSec

rvMeasurement.consistencyRangeBiasM
rvMeasurement.consistencyVelocityBiasMps
rvMeasurement.matchRangeResidualM
rvMeasurement.matchVelocityResidualMps
rvMeasurement.matchRangeThresholdM
rvMeasurement.matchVelocityThresholdMps
rvMeasurement.matchNormalizedScore

rvMeasurement.preSpeedGateCandidateCount
rvMeasurement.postSpeedGateCandidateCount
rvMeasurement.postOneToOneCandidateCount
```

这些字段可以用来检查：

- 哪些峰的有效 SINR 较低；
- 不同目标实际使用了多宽的 R/V 门；
- 候选被拒绝是因为距离残差还是速度残差；
- 速度硬门和一对一约束分别剔除了多少候选；
- 当前结果是否由跨频段固定偏差主导。

## 16. 参数调整建议

### 16.1 `sigmaMultiplier`

- 减小：降低假警，但增加漏检风险；
- 增大：降低漏检，但允许更多组合候选；
- 建议在包含不同 SNR、目标间距和漏峰情况的 Monte Carlo 数据上确定。

### 16.2 `noiseEstimator`

- 多目标、强旁瓣场景优先使用 `localCfar`；
- 研究纯热噪声 CRLB 时可使用 `robustOffTarget`；
- 两种噪声估计都会保存到 `rvMeasurement`，便于对照。

### 16.3 `delayStdFloorBins`

- 该值描述 CRLB 之外的峰拟合和模型精度下限；
- 过小会在高 SNR 下形成不现实的窄门；
- 过大会使门限失去随 SNR 收紧的能力。

### 16.4 `maxPeakRankDifference`

- 峰数一致且目标距离不交叉时使用 `0`；
- 可能漏峰时使用 `1` 或更大值；
- 关闭排序约束时必须依赖更强的跨帧跟踪或全局数据关联，避免组合假目标重新出现。

## 17. 适用范围和限制

当前推导采用以下假设：

1. 三个频段的随机时延误差近似独立；
2. 峰附近误差可用局部高斯近似描述；
3. LFM 频谱近似矩形；
4. `rankMedian` 假设多数频段峰具有相同的距离排序；
5. 单个 CPI 内目标距离顺序没有发生交叉；
6. 当前门限不代替跨 CPI 跟踪、航迹关联或运动先验。

如果后续加入真实噪声、强杂波、不同频段漏检或目标距离交叉，应进一步考虑：

- 使用完整的 \((\Delta R,\Delta V)\) 协方差矩阵和马氏距离门；
- 对不同频段时延误差建立相关项；
- 使用全局最优的一对一分配代替贪心分配；
- 跨 CPI 使用 Kalman/IMM 等运动预测门；
- 通过 Monte Carlo 仿真标定 CRLB 与三点峰拟合实际方差之间的比例。

## 18. 代码位置

- 参数入口：`code/multiuav_iasr_fullpath.m`
- SINR、CRLB、雅可比传播和候选关联：`code/multiuav_rvpath.m`
- 双频段 R/V 解算：`radar.Cal_rv_mat`
- legacy 固定门参考：`code/time_mul_tar_Meas_6LFM_opt.m`

