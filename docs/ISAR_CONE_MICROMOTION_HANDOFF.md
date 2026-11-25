# 锥体微动 ISAR 分析交接

更新：2026-10-21。本文对应本次锥体微动分支；基于 RadarProject 原提交 `33da190`，不包含工作区内其他单脉冲、卫星或硬件实验的未提交改动。

## 1. 从哪里开始

当前主入口是 `experiments/run_isar_tip_motion_retest.m`：真实公共库宽带回波、每 CPI 固定本振、原 ISARpath 前端、匀加速先验复相位校正、锥尖距离门、慢时间频谱切片。最新交付是 **14 个物理场景、每场景 50 个 CPI、0–3 s 显示**。

其他入口按研究阶段保留，不应把不同阶段的补偿链和结果混用：

| 任务 | 入口 |
|---|---|
| 原始 fullpath 接口，切换自旋/进动与相位分支 | `code/multiuav_iasr_fullpath.m`、`code/multiuav_isarpath.m` |
| 四点 32 ms 方位功率 −3 dB 全宽实验 | `run_isar_spin_precession_32ms` |
| 静止目标 5 s 四种微动、1/4 kHz 对照 | `run_isar_long_micro_motion` |
| 径向运动、匀速先验的早期解析 IF 实验 | `run_isar_moving_cv_micro_motion` |
| 实际 lib RF 与固定本振独立混频的 CV/CA 对照 | `run_isar_lib_fixed_lo_study` |
| 32 个物理工况的完整参数扫描 | `run_isar_motion_grid_suite` |
| 原 ISARpath 全链、两种典型场景、整数/分数距离对齐 | `run_isarpath_two_scene_sliding` |
| 1/3 Hz 小锥角、三种补偿方式、30 号字动图 | `run_isarpath_small_cone_font30` |
| 最新 14 场景的锥尖慢时间分析 | `run_isar_tip_motion_retest` |
| 仅重绘 14 张融合颜色图 | `render_isar_tip_merged_colours` |
| ±3°/±5°摆动可见性与 1 Hz 分量分析 | `analyze_isar_tip_swing_visibility` |
| 锥体结构/LOS 示意图 | `render_isar_cone_structure_font30`、`render_isar_cone_geometry_los` |
| 自旋及自旋加进动的 1 s 模型动图 | `render_isar_micro_motion_gifs` |

表中的 MATLAB 函数位于 `experiments/`。归因诊断、频率审计和汇总脚本也在该目录；对应已有测试位于 `tests/test_isar*.m`、`tests/verify_isar*.m`。历史逐点成像说明另见 [2026-10-15 实验汇总](ISAR_SPIN_PRECESSION_EXPERIMENT_SUMMARY_20261015.md)。

## 2. 依赖与启动

- 本机验证环境：Windows、MATLAB R2025b；Signal Processing Toolbox 提供 `resample`。`NativeWorkers > 1` 需要 Parallel Computing Toolbox，必须使用进程池，不要改成线程池。
- 公共库仓库：<https://github.com/ragiko0104/matlab-research-lib>。
- 本次核对的公共库提交：`c80bbe21af9a579fddc62cc1e56258c0f483a041`。其工作区在交接时干净，相关依赖已被该仓库跟踪。
- 关键公共接口：`scatterermodel.createLegacyDiskTarget`、`scatterermodel.createConicalMicroMotion`、`scatterermodel.evaluateTarget`、`radar.echo_generate_target_lfm`、`radar.lfm_local_reference_beat`、`radar.doppler_dispersion_cmp`、`radar.phase_cmp`、`radar.ISAR_imaging` 以及 `signal` 包。
- 本分支不复制公共库，不依赖未提交的项目 `lib/+scatterermodel/createStarlinkV2MiniTarget.m`。
- 中文图默认 `Microsoft YaHei`。在其他平台需要安装可用中文字体，或明确修改绘图字体；字号要求仍是 30。

推荐保持相邻布局，因为多个历史实验入口直接加入相邻公共库路径：

```text
MatlabProj/
├── RadarProject/
│   ├── code/  functions/  experiments/  tests/  docs/
│   └── results/                  # 运行生成，被 .gitignore 排除
└── matlab_research_lib/
    └── lib/
```

在 RadarProject 根目录启动 MATLAB：

```matlab
root = pwd;
libRoot = fullfile(fileparts(root),'matlab_research_lib');
assert(isfolder(fullfile(libRoot,'lib')));
addpath(fullfile(root,'functions'),fullfile(root,'experiments'), ...
    fullfile(root,'tests'),fullfile(libRoot,'lib'));
assert(~isempty(which('radar.echo_generate_target_lfm')));
assert(~isempty(which('scatterermodel.createConicalMicroMotion')));
maxNumCompThreads(2);
```

`startup_project.m` 支持 `MATLAB_RESEARCH_LIB_ROOT`，但并非所有历史入口都统一解析此变量。若使用非相邻布局，先显式加入正确的公共库 `lib/`；不要只设置环境变量就假设所有入口都完成了初始化。

## 3. 物理模型与口径

四个理想点散射体：P1 为锥尖，P2–P4 为底部三点。尖端到基底轴向间距 1.5 m，底面半径 0.3 m，相干系数 `[1;0.3;0.3;0.3]`。点散射体没有有限物理直径。旋转中心是原模型体坐标原点；P1 坐标为 `(0,-0.75,0)`，不是固定在世界坐标中的进动支点。

初始姿态为世界 Z 轴转 −135°；雷达至目标 LOS 为世界 +X。无额外摆动时，固定进动轴/纯自旋轴与 LOS 的锐角为 45°。

```text
纯自旋：       C(t) = C0 Ry(2π fs t)
自旋加进动：   C(t) = C0 Ry(2π fp t) Rz(alpha) Ry(2π fs t)
叠加整体摆动： 在上述姿态左侧乘 Rz(beta(t))
beta(t) = -A cos(2π fw t)，对称范围 [-A,+A]，初始角在最小值
```

这里约定的“章动”是**整个锥体绕世界固定 Z 轴往复摆动**，不是进动半锥角随时间变化。`fs` 是最后一个欧拉旋转因子的频率；无摆动时，物理轴向角速度投影对应 `fs + fp*cos(alpha)`。P1 在自旋轴上，纯自旋不改变其相对质心位置；因此不能从 P1 轨迹辨识自旋频率，也不能把 P1 运动频率直接写成 `fs+fp`。

方位定标采用目标/窗口各自的局部 LOS 有效转动率，见 `isar_cross_range_geometry.m`；进动/摆动时成像平面会变化，不能把局部等效距离轴当成全程固定的精确物理横向坐标。点响应指标统一用**功率 −3 dB 全宽**，不使用主瓣零点宽度或 90% 能量宽度；零填充不增加真实分辨率。

## 4. 最新 14 场景与处理链

公共参数：5 s 观测、200 MHz 去斜后采样、11.8–17.8 GHz LFM、100 μs 脉宽、载频 14.8 GHz、`K=6e13 Hz/s`、每 100 ms 一个窗口、PRF 1 kHz、32 脉冲/名义 32 ms CPI。首末脉冲起点跨度为 31 ms。目标中心 `R(t)=10050+3400t+50t² m`；加速度为代码取值 **100 m/s²**，不要把它写成精确的标准 10 g。

| 组 | 自旋 fs/Hz | 进动 fp/Hz | 半锥角 alpha | 固定 Z 轴摆动 | 场景数 |
|---|---|---|---|---|---:|
| 1 | 3 | 3 | 1°、3° | 无 | 2 |
| 2 | 2、3.5 | 2、3.5 | 5°、7° | 无；全部交叉组合 | 8 |
| 3 | 2 | 3.5 | 5°、7° | 1 Hz，±3°、±5°；全部交叉组合 | 4 |

最新流程如下：

1. `isarpath_synthesize_profiles` 在 1 GHz 内部快时间网格调用公共库生成实际解析宽带目标回波；与每 CPI 固定本振共轭混频，经过原库的正频支路、滤波、重采样和快时间去直流，得到 200 MHz IF。内部网格是数值合成约定，不代表真实 RF ADC 能直接采样 6 GHz 带宽。
2. `isarpath_extract_segment`、`isarpath_range_profiles` 复用原 ISARpath 时域截取、色散补偿、频域移窗/截取。**最新 CA 分支仍共用原 CV 包络预测、色散和截窗链**，再在原相位校正位置使用 CA 运动先验。
3. 先用四点相干场景求共同包络互相关距离对齐，再把相同位移用于诊断通道。最新锥尖实验用分数移位，无特显点补偿，无慢时间去均值。
4. `isarpath_prior_phase` 只使用外部给定的距离、速度、加速度及处理坐标，计算共同复相位因子；不使用真实散射点微动或实际场景峰值相位。`MODEL_CV` 不使用加速度。它不仅移动频域窗，也计算正 IF 相位、色散相位及截窗/对齐坐标的相位参考。
5. 距离门从**完整四点场景**选择。搜索区相对于先验预测的处理后中心为 `[-0.70,-0.36] m`，门半宽 3 个距离 bin。不可改回相对图像坐标零点搜索，否则晚期 CPI 会选错门。孤立 P1 通道只在选门后审计主峰覆盖及尾点泄漏，不参与选门或补偿。
6. 对门内 HRRP 各 bin 做慢时间 FFT，再非相干累加功率并插值提取峰值。64 倍 FFT 零填充仅用于频率插值。保存小型锥尖频谱切片/轨迹/审计表，不保存 RF、IF 或完整 HRRP。

不要添加独立的“移位相位恢复”因子。对称域分数循环移位本身的性质见 `isar_symmetric_fractional_shift.m` 和相应测试；先验相位函数已经按实际截窗/对齐坐标定义参考。历史关于“整数移位必然引入额外相位台阶”的宽泛解释不能直接作为当前流程的结论。

## 5. 运行、恢复与绘图

完整运行最新 14 场景：

```matlab
run_isar_tip_motion_retest(struct( ...
    'Groups',[1 2 3], ...
    'ThirdConingHalfAnglesDeg',[5 7], ...
    'ThirdSwingAmplitudesDeg',[3 5], ...
    'NativeWorkers',8));
```

`Groups` 默认只含 1、2，第三组必须显式传入两个角度集合。首组完成后先检查 100 ms 窗口更新的采样条件，通过后才进入后组。本次 P1 主要位置运动频率为 3 Hz，小于 10 Hz 切片采样率的 5 Hz 奈奎斯特上限；首组的峰值瞬时多普勒另为约 17.215/51.623 Hz，与“运动调制频率”不是同一个量。

默认输出 `results/isar_tip_motion_retest_20261020/`，包含三个 `groupN/` 子目录。每场景 `<caseId>_tip_result.mat` 是可恢复的小型检查点，主要保存 1025×50 的锥尖门频谱功率、轨迹和几何/波形审计。场景阶段结束后写 `scene_similarity.csv`；全部完成后写 `all_scene_similarity.csv`、`all_scene_figures.csv`、`batch_summary.mat`。

完整 native 运行会比较耗时。本机 8 个进程工作者的该次运行约半小时；这只是历史实测量级，不是其他机器的时延保证。可用 `NativeWorkers=1` 串行运行。只生成被选 CPI 内的脉冲，不生成不参与滑窗的间隔样本。

检查点基于固定场景 ID 等有限字段恢复，不是任意参数修改后的可靠缓存失效系统。改变模型/参数时给 `OutputDirectory` 指定新目录；不要沿用旧目录制造伪复测。单组未完成时，组级 CSV 可能尚不存在，应查看逐场景检查点。

仅重绘（必须先有这些小型结果；Git 不包含 MAT/PNG/FIG）：

```matlab
render_isar_tip_merged_colours();
analyze_isar_tip_swing_visibility();
```

最终 14 张颜色图使用原文件名 `<caseId>_colour_3s.png/.fig`，覆盖原颜色图版本：30 号字、0–3 s、色标 −10～0 dB、背景不隐藏、无轴标题；白虚线为物理理论，青实线为原始提取峰值的谐波拟合，圆点为各 CPI 峰值。标题第二行或第三行显示原始峰值相对于物理理论的 `R`。下方图例透明、文字黑色；白色理论曲线在图例中使用深灰虚线样本以保证可见。

拟合用全部 50 个原始峰值，基函数为常数、线性项以及 0.5:0.5:4.5 Hz 的正余弦，共 20 参数。理论轨迹不参与拟合；R 仍由前 3 s 原始峰值计算。该固定频率网格适用于本组实验，不是任意未知频率下的通用运动识别算法，也不能恢复采样丢失的信息。现有单独理论对比图保持独立保留。

`analyze_isar_tip_swing_visibility` 使用同参数的无摆动**仿真**轨迹作基线，统一纵轴显示摆动增量，输出 `tip_swing_visibility_3s.png/.fig`。这是诊断对照，不是实际单段观测中可直接获得的无摆动基线。

## 6. 已有结果与解释边界

本分支附带 `docs/isar_micro_motion_reference/` 下四份小型 CSV：14 场景相似度、首组采样门控、摆动分量幅度、摆幅对比。它们是本地完整运行的归档数值，不表示本次发布又重新运行了全部 700 个 CPI。它们不含原始回波，也不含机器绝对路径。完整图像和 MAT 仍只在本地 `results/`，需要重跑才能在新克隆中生成。

- 本次 14 场景共 700 CPI；前 3 s 各用 30 点计算 R，范围约 **0.999939–0.999991**。所有 CPI 审计均覆盖 P1 主峰；最大原始拍频约 **69.439 MHz**，小于 100 MHz，四点相干和误差最大约 `6.18e-16`。
- 原始峰值相对纯物理多普勒仍有约 5 Hz 偏移，前 3 s 各场景平均约 4.967–5.044 Hz。共同质心相位校正后还保留逐点 RVP/距离耦合项，近似 `-4*K*v_center*DeltaR/c²`。历史 5 s 平均约 5.2 Hz；不要把它解释为运动混叠频率，也不要为了提高 R 偷做时间平移/偏置扣除。
- R 衡量各自处理与各自理论的吻合，不是不同摆幅的可分性，更不能证明 5 cm 空间分辨率。Pearson 相关系数对整体幅度缩放和常量偏移不敏感。
- 在相同进动半锥角下，±3°/±5°摆动的 1 Hz 提取分量分别约 17.09/28.47 Hz（半锥角 5°）和 17.03/28.37 Hz（半锥角 7°）。进动 3.5 Hz 分量仍约 99/138 Hz，前 3 s 两摆幅轨迹最大差别仅约 12.51/13.09 Hz。原图自动纵轴和较宽谱条带会弱化差异；32 ms 矩形窗匀频单点的理想频谱 −3 dB 全宽约 `0.885893/0.032 = 27.684 Hz`。
- PRF 对应的慢时间多普勒奈奎斯特范围与切片更新频率必须分别审计。1 kHz PRF 下尾点可能混叠；距离门显著压低其对 P1 的泄漏，不代表完整四点图已经无混叠。4 kHz 对照采用 128 脉冲保持名义 32 ms。
- 第三组含 1 Hz 摆动及其与进动的边带；远场解析核查存在很弱的 5.5 Hz 高阶项，100 ms 切片不能严格无混叠地保留所有高次边带。首组通过 3 Hz 检查不等于所有非线性频率分量都严格带限。
- 一维共同包络对齐不能同时消除四点各自的差异距离徙动或非线性微动相位。按真实逐点轨迹的 oracle 校正仅用于归因，不是实际可用的补偿分支。
- 当前实验无噪声、无先验误差；±3°/±5°识别结果不能直接推广成未知几何、未知噪声条件下的摆角估计精度。

## 7. 验证与后续工作

本次发布的实际检查、通过项和未重跑范围见 [发布验证记录](ISAR_CONE_MICROMOTION_VALIDATION.md)。可以独立运行以下现有测试；先按第 2 节加入路径：

```matlab
test_isar_micro_motion;
test_isar_fractional_alignment;
test_isarpath_full_chain;
test_isarpath_native_frontend;
test_isarpath_script_integration;
```

`test_isarpath_prior_phase_dense`、`test_isar_motion_grid` 等较重测试可按需运行。部分旧 `test_isar_*` 和 `verify_isar_*_outputs` 明确读取先前实验结果，应先运行对应生产入口；不应把缺少被 Git 忽略的 MAT 文件误认为新克隆的算法错误。

建议继续：固定相同纵轴对比摆动参数；从 1 Hz 分量和边带进行含 LOS 几何的参数拟合，并加入噪声/先验误差；若要辨识自旋，纳入尾点并先处理 PRF 混叠；若要证明 5 cm 可分辨，还需相邻双点验证，不能只看单点图或 R。

本次 Git 只发布微动相关源码、已有测试、交接和小型参考 CSV；不发布 `results/`、缓存、MAT/FIG/PNG、原始 RF/IF/HRRP、压缩包或其他工作区实验。
