# RV_ISAR_Path

本分支（`RV_ISAR_Path`）为从原 `RadarProject` `main` 分支清洗迁移而来的**测距测速（RV）+ ISAR 成像核心仿真链路**，已移除全部实验室硬件/FPGA/XSim/AWG/TCP/上位机脚本与旧版重复本地函数，所有底层信号与雷达算子统一依赖公共库 [`matlab-research-lib`](https://github.com/ragiko0104/matlab-research-lib)。

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
│   ├── multiuav_rvpath.m                          # 六段 LFM 测距测速链路（全 lib 驱动）
│   └── multiuav_isarpath.m                        # 宽带去斜 ISAR 成像与质量评估链路（全 lib 驱动）
├── functions/
│   └── run_uav_existing_pipeline.m                # 透传至 uavmodel.scenarios.runExistingRadarPipeline
├── docs/
│   └── fullpath_rv_gate_formula.md                # 六段 LFM SNR/SINR 一致性门限推导与配置说明
└── tests/
    ├── test_rvcal_fpga_alignment.m                # CFAR/Top-K 阈值与零基时延回归测试
    ├── test_legacy_four_scatterer_unification.m   # 四点模型几何与 LOS 投影一致性测试
    └── test_rv_isar_fullpath_smoke.m              # 10 目标 RV + ISAR 全链路无人值守回归测试
```

## 快速开始

### 1. 运行 RV + ISAR 全链路（默认四点或 UAV 模型）

```matlab
startup_project
run_RVpath = true;
run_ISARpath = true;
targetModelTypes = repmat("legacyDisk", 10, 1); % 可选 "uav" / "legacyDisk" / "rigidPointCloud" / "point"
run("code/multiuav_iasr_fullpath.m")
disp(isarQuality.summary)
```

### 2. 运行紧凑 UAV 场景 A–E

```matlab
targetModel = "uavScatterer";
uavScenario = "ALL";       % A, B, C, D, E, or ALL
uavEchoMode = "stopAndGo"; % 或 "fastTimeExact"
main
```

### 3. 运行回归测试

```matlab
run("tests/test_rvcal_fpga_alignment.m")
run("tests/test_legacy_four_scatterer_unification.m")
run("tests/test_rv_isar_fullpath_smoke.m")
```
