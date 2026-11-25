# RV_ISAR_Path

本分支（`RV_ISAR_Path`）为从原 `RadarProject` 清洗并升级而来的**测距测速（RV）+ 模块化 ISAR 成像基础链路**，作为 `Radar_Repo` 中共享的 RV/ISAR 唯一代码源（Single Source of Truth）：

- 已将**仿真参数、目标个数、目标运动参数、目标散射体模型、波形参数与各类开关**全部抽离到独立的 `config/rv_isar_default_config.m`，在 `fullpath` 外部通过 `simCfg.switches.runMode` 即可灵活控制只跑 `"rv"`、只跑 `"isar"` 或 `"both"`；
- 已移除全部硬件/FPGA/XSim/AWG/TCP/上位机脚本与旧版重复本地函数，底层算子统一依赖 [`matlab-research-lib`](https://github.com/ragiko0104/matlab-research-lib)；
- 已同步合入向后兼容的模块化 `isarpath_*` 算子层，供 `MD_ZT` 分支直接复用。

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
├── main.m                                         # 外部参数配置与运行总入口
├── config/
│   └── rv_isar_default_config.m                   # 独立参数/目标/散射体/波形/开关配置模板
├── code/
│   ├── multiuav_iasr_fullpath.m                   # 纯执行调度脚本（读取外部 simCfg）
│   ├── multiuav_rvpath.m                          # 六段 LFM 测距测速链路（固定 fastTimeExact，全 lib 驱动）
│   └── multiuav_isarpath.m                        # 模块化宽带去斜 ISAR 成像与质量评估链路（全 lib 驱动）
├── functions/
│   ├── isarpath_*.m                               # 7 个共享 ISARpath 模块化算子
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
    └── test_rv_isar_fullpath_smoke.m              # 外部 simCfg 三种 runMode("rv"/"isar"/"both") 回归测试
```

## 快速开始

### 方式 1：直接编辑并运行 `main.m`
在 `main.m` 中修改 `simCfg.switches.runMode`（`"rv"` / `"isar"` / `"both"`）、目标参数、散射体模型或波形参数后直接运行 `main`。

### 方式 2：在任意外部脚本中构造 `simCfg` 并调用 `fullpath`

```matlab
startup_project
simCfg = rv_isar_default_config();

% 选择仅跑 RV_path、仅跑 ISAR_path，或两者串联
simCfg.switches.runMode = "both";          % "rv" | "isar" | "both"

% 配置目标个数、运动参数与散射体模型
simCfg.targets.targetIndices = [1 2];      % 或直接指定 simCfg.targets.R0 / v_Mach / a_g
simCfg.scatterer.modelTypes  = "legacyDisk";

% 配置波形与成像选项（仅 ISAR 支持 stopAndGo）
simCfg.isar.Na       = 32;
simCfg.isar.echoMode = "fastTimeExact";

run("code/multiuav_iasr_fullpath.m")
```

### 运行回归测试

```matlab
run("tests/test_rvcal_fpga_alignment.m")
run("tests/test_legacy_four_scatterer_unification.m")
run("tests/test_rv_isar_fullpath_smoke.m")
```
