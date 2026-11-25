# 锥体微动分支发布验证

日期：2026-10-21。说明与入口见 [handoff](ISAR_CONE_MICROMOTION_HANDOFF.md)。

## 环境与范围

- Windows、MATLAB R2025b；公共库固定为 `matlab-research-lib` 提交 `c80bbe21af9a579fddc62cc1e56258c0f483a041`。
- 测试在独立 Git worktree 中运行，基于 RadarProject `33da190` 加本次明确复制的微动源码；没有使用原工作区其他单脉冲/卫星/硬件实验的未提交改动。
- 在运行测试前断言 `isar_apply_micro_motion`、`isarpath_prior_phase` 来自待发布检出，`radar.echo_generate_target_lfm` 来自上述公共库。
- 该检出位于本机 `.tmp/`，不满足历史入口的默认相邻目录布局；已显式加入真实公共库路径。历史入口因此产生过无效相邻路径警告，但函数解析检查及以下测试均通过。按 handoff 推荐的相邻仓库布局可避免这些警告。

## 本次实际运行

| 检查 | 结果 | 耗时约/s | 验证内容 |
|---|---|---:|---|
| `test_isar_micro_motion` | PASS | 0.552 | 四种运动的向量化状态与通用库状态一致、质心运动保持、P1 纯自旋静止、局部方位定标 |
| `test_isar_fractional_alignment` | PASS | 0.122 | 对称域 Fourier 位移、整数极限、可逆性、分数互相关对齐和有限孔径点响应 |
| `test_isarpath_full_chain` | PASS | 0.220 | 默认链与冻结旧表达式一致、距离对齐、特显点及先验相位分支位置 |
| `test_isarpath_native_frontend` | PASS | 7.059 | 实际库 RF 的缓存/非缓存一致性、四点相干和、独立 CA 中心参考、CV 不泄漏加速度 |
| `test_isarpath_script_integration` | PASS | 6.428 | 实际 `multiuav_isarpath.m` 与流式全链实验的 HRRP/相位后结果/图像一致 |
| `test_isar_cv_micro_motion` | PASS | 1.643 | 移动目标 RF/IF 参考、色散补偿、未知加速度下的 CV 对照 |
| `test_isar_long_micro_motion` | PASS | 1.033 | 长观测解析 IF 与独立 RF 乘积、快时间状态、纯自旋锥尖和共同 RCMC 回放 |
| 最新锥尖流程单 CPI 验证 | PASS | 7.3 | 真实 lib 生成 32 脉冲、1 kHz、32 ms CPI，分数对齐、CA 先验、锥尖距离门与输出约束 |

实际父脚本与流式前端的 HRRP 相对误差约 `1.74e-16`。native 前端测试还验证了中心点 CA 校正后残余相位小于 0.02 rad、CV 分支不会因输入加速度字段变化而改变相位因子。

最后一项使用自旋 3 Hz、进动 3 Hz、半锥角 1°的一个 CPI，`NativeWorkers=1`。保存的 P1 主峰在选定距离门内；断言最大拍频低于 100 MHz、四点相干和误差低于 `1e-11`，并断言没有保存原始回波/完整 HRRP、没有开启特显点补偿或慢时间去均值。它验证生产入口闭环，不替代完整 5 s 性能验证。

## 本次未重跑的范围

- 没有重新生成全部 14 场景的 700 CPI，也没有重新执行历史 32 工况全网格或所有逐点归因控制实验。
- `docs/isar_micro_motion_reference/` 的四个 CSV 来自此前完成的本机完整运行；handoff 中的相关系数、拍频上限和摆动幅度结论据此归档。
- 部分旧测试要求先前结果 MAT。没有把缺少未提交的结果文件列为测试失败，也没有声称所有 `tests/` 文件均已运行。
- 本次发布未改变运动模型、RF 合成、补偿或成像参数。对 `run_isar_motion_grid_suite` 增加了首次运行记录 CSV 不存在时的创建保护；未重新运行完整 32 工况来触发最后的记录写入阶段。
- 完整 RF/IF/HRRP、结果 MAT/FIG/PNG、缓存和本机临时验证脚本均不纳入 Git。

提交前另检查明确文件清单及 `git diff --cached --check`；发布后核对远端分支与本地提交一致。
