# RV_ISAR_Path

本分支（`RV_ISAR_Path`）为从原 `RadarProject` 清洗并升级而来的**测距测速（RV）+ 模块化 ISAR 成像基础链路**，作为 `Radar_Repo` 中共享的 RV/ISAR 唯一代码源（Single Source of Truth）：

- 已移除全部硬件/FPGA/XSim/AWG/TCP/上位机脚本与旧版重复本地函数；
- 底层信号与雷达算子统一依赖公共库 [`matlab-research-lib`](https://github.com/ragiko0104/matlab-research-lib)；
- 已同步合入锥体微动研究中沉淀的**向后兼容模块化 `isarpath_*` 算子层**（支持分数距离对齐 `FRACTIONAL`、`MODEL_CV`/`MODEL_CA` 先验相位补偿及 `ISAR_MOTION_OPTIONS` 姿态接入），使 `MD_ZT` 分支无需复制任何 RV/ISAR 基础脚本。

## 依赖配置（`matlab_research_lib`）

`startup_project.m` 按以下优先级自动定位 `matlab_research_lib/lib`：
1. `local_config.m` 或当前工作区变量 `MATLAB_RESEARCH_LIB_ROOT`
2. 环境变量 `MATLAB_RESEARCH_LIB_ROOT`
3. 工程内 `external/matlab_research_lib/lib`
4. 同级目录 `../matlab_research_lib/lib`
5. 默认本机路径 `D:\MatlabProj\matlab_research_lib\lib`

如需自定义路径，复制 `local_config_template.m` 为 `local_config.m` 并修改 `MATLAB_RESEARCH_LIB_ROOT` 即可。

## 目录结构

```text
RV_ISAR_Path/
├── startup_project.m                              # 路径初始化与公共库定位
├── local_config_template.m                        # 本机路径覆盖模板
├── main.m                                         # 总入口（fullpath / uavScatterer）
├── code/
│   ├── multiuav_iasr_fullpath.m                   # 多模型场景配置与 RV + ISAR 串联入口
│   ├── multiuav_rvpath.m                          # 六段 LFM 测距测速链路（固定 fastTimeExact，全 lib 驱动）
│   └── multiuav_isarpath.m                        # 模块化宽带去斜 ISAR 成像与质量评估链路（全 lib 驱动）
├── functions/
│   ├── isarpath_settings.m                        # 共享 ISARpath 后处理配置契约
│   ├── isarpath_extract_segment.m                 # 接收时域截窗
│   ├── isarpath_range_profiles.m                  # 色散补偿 + 距离 FFT + 移窗截取
│   ├── isarpath_align_profiles.m                  # 整数/分数包络互相关距离对齐
│   ├── isarpath_prior_phase.m                     # MODEL_CV / MODEL_CA 解析先验复相位参考
│   ├── isarpath_phase_image.m                     # 相位补偿与 radar.ISAR_imaging 成像
│   ├── isarpath_synthesize_profiles.m             # 流式逐 CPI 原生宽带回波与去斜前端
│   ├── isar_physical_fractional_shift.m           # 物理原点分数循环移位
│   ├── isar_apply_micro_motion.m                  # 微动姿态接入（封装 scatterermodel.createConicalMicroMotion）
│   ├── isar_rigid_micro_motion_state.m            # 刚体微动矢量化状态解算
│   ├── isar_cross_range_geometry.m                # 局部 LOS 有效转动率与横向定标
│   └── run_uav_existing_pipeline.m                # 透传至 uavmodel.scenarios.runExistingRadarPipeline
├── docs/
│   └── fullpath_rv_gate_formula.md                # 六段 LFM SNR/SINR 一致性门限推导与配置说明
└── tests/
    ├── test_rvcal_fpga_alignment.m                # CFAR/Top-K 阈值与零基时延回归测试
    ├── test_legacy_four_scatterer_unification.m   # 四点模型几何与 LOS 投影一致性测试
    └── test_rv_isar_fullpath_smoke.m              # RV (fastTimeExact) + ISAR (stopAndGo) 串联回归测试
```

## 快速开始

```matlab
startup_project
run_RVpath = true;
run_ISARpath = true;
targetModelTypes = repmat("legacyDisk", 10, 1);
run("code/multiuav_iasr_fullpath.m")
disp(isarQuality.summary)
```

### 运行回归测试

```matlab
run("tests/test_rvcal_fpga_alignment.m")
run("tests/test_legacy_four_scatterer_unification.m")
run("tests/test_rv_isar_fullpath_smoke.m")
```
