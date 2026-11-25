# MD_ZT（锥体微动测量与识别）

本分支（`MD_ZT`）为从原 `RadarProject` 的 `codex/cone-micro-motion-analysis-20261021` 分支清洗迁移而来的**锥体微动测量、慢时间轨迹提取、盲模式识别（H0–H3）、尾点自旋识别及多帧拟合加速度 CA 成像链路**。

## 清洗与分层架构

- **已整块移除**：原分支误带入的 131 个三态解混 / 交错 UAV / FDM / MZM 光子前端 / 单脉冲编码文件、35+ 个阶段性排障/早期固定本振/PS1 汇总/论文渲染与 PPT 脚本，以及 `references/fpga_golden_0ef90c0/` 下的全部 FPGA 脚本（仅保留纯参数配置 `config/golden_10targets_config.json`）。
- **零重复分层依赖**：
  1. 底层雷达与信号算子依赖 [`matlab-research-lib`](https://github.com/ragiko0104/matlab-research-lib)；
  2. 共享 RV + ISAR 主链路（`multiuav_iasr_fullpath.m`、`multiuav_rvpath.m`、`multiuav_isarpath.m`、`isarpath_*.m`、`isar_apply_micro_motion.m` 等）直接复用同级目录的 **`RV_ISAR_Path`** 分支，本分支不保留重复副本。

## 依赖配置

`startup_project.m` 自动加载：
1. `matlab_research_lib/lib`（优先 `local_config.m` / `MATLAB_RESEARCH_LIB_ROOT` → `external/matlab_research_lib` → `../matlab_research_lib` → `D:\MatlabProj\matlab_research_lib`）
2. `RV_ISAR_Path`（优先 `local_config.m` / `RV_ISAR_PATH_ROOT` → `../RV_ISAR_Path` → `D:\MatlabRepo\RV_ISAR_Path`）

如需自定义路径，复制 `local_config_template.m` 为 `local_config.m` 修改即可。

## 核心入口一览

| 任务 | 入口文件 |
|---|---|
| 锥体微多普勒 + ISAR 基础仿真（全 `lib` 驱动） | `code/DT_MD_simulation.m` |
| 14 场景原生宽带回波锥尖慢时间提取 | `experiments/run_isar_tip_motion_retest.m` |
| 锥尖轨迹盲微动模式识别（H0–H3 嵌套 F 阶梯） | `experiments/identify_isar_tip_micro_motion.m` |
| 尾部散射点四轨联合自旋识别 | `experiments/identify_isar_spin_from_tail_points.m` |
| 全场多帧测距拟合加速度 CA → 锥尖识别端到端链路 | `experiments/isar_chain_identify_ch1.m` |
| 10 目标多帧拟合加速度 CA 成像 | `experiments/isar_multi_frame_10targets_v2.m` |
| 锥尖二维分辨率计量对照（含无微动理想点对照） | `experiments/measure_tip_2d_resolution.m` |
| 测距噪声加速度拟合验证 | `experiments/isar_range_fit_noise.m` |
| 灵敏度 / 带宽 / 积分时间 / 噪声与分辨率扫描 | `experiments/sweep_isar_*.m`、`test_isar_integration_*.m`、`test_isar_micro_motion_resolution.m`、`identify_isar_tip_noise_study.m` |

## 回归测试

```matlab
startup_project

% 共享 ISARpath 与跨仓库对接回归
test_isar_micro_motion
test_isar_fractional_alignment
test_isarpath_full_chain
test_isarpath_native_frontend
test_isarpath_script_integration

% v14 盲模式识别与尾点自旋回归
test_isar_cone_point_kinematics
test_isar_tip_motion_mode_identification
test_isar_spin_from_tail_points

% PR #2 拟合加速度 CA 契约与理想点台架回归
test_isar_cone_point_input_validation
test_isar_chain_prior_isolation
test_isar_golden_lo_offsets
test_isar_tip_to_coherent_ratio
test_isar_ca_run_contract
test_ideal_point_bench_ca
```
