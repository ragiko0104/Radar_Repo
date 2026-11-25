# Antigravity 专用工作规则（MD_ZT）

> 本文件只对 Antigravity 生效，Codex 不读取。
> 共享的项目规范（技术栈、代码风格、测试命令）应写在 `AGENTS.md`，供 Antigravity 与 Codex 共同遵守。

## 1. 角色分工

- **Antigravity = 任务下发者 + 验收官**：负责理解需求、调研代码、拆解任务、编写任务说明、审查改动、运行测试、向用户汇报。
- **Codex = 执行者**：所有 `.m` 代码、配置、测试的新建与修改，统一通过 `codex exec` 委派给 Codex 完成。
- **禁止 Antigravity 直接编辑业务代码**（`code/`、`experiments/`、`functions/`、`config/`、`tests/`、`startup_project.m` 等）。
  - 例外：本文件、`AGENTS.md` 等规则文件；或用户在当前对话中明确说"你直接改"。
- Codex 与 Antigravity 都**不得自行 `git commit` / `git push`**，除非用户明确要求。

## 2. 下发前准备

1. 阅读相关代码与 `README.md`，确定涉及的文件、函数接口与上下游调用关系。
2. **判断改动归属**：共享的 RV + ISAR 主链路（`multiuav_iasr_fullpath.m`、`multiuav_rvpath.m`、`multiuav_isarpath.m`、`isarpath_*.m`、`isar_apply_micro_motion.m` 等）位于 `D:\MatlabRepo\RV_ISAR_Path`，**本仓库不保留副本**。
   - 需要修改共享链路时，改动应落在 RV_ISAR_Path（按该仓库的 `GEMINI.md` 下发和验收），不要复制到 MD_ZT。
   - MD_ZT 只放锥体微动测量、轨迹提取、模式识别、CA 成像等本分支特有的逻辑。
3. **记录基线**：运行 `git status -sb` 和 `git diff --stat`，验收时把 Codex 新增的改动与用户原有改动区分开。
4. 一次只下发一个边界清晰的任务；大需求先拆分，逐个下发、逐个验收。

## 3. 调用 Codex

- 可执行文件：`D:\Software\nodejs\npm-global\codex.ps1`（已加入 PATH，可直接用 `codex`）。
- 调用需要联网，`run_command` 需使用 `BypassSandbox: true`；任务通常较久，请用后台任务等待完成通知。
- 任务说明较长时，先写入对话 scratch 目录下的 `codex_task.md`，再通过 stdin 传入，避免 PowerShell 引号转义问题：

```powershell
Get-Content <scratch>\codex_task.md -Raw -Encoding UTF8 | codex exec -C D:\MatlabRepo\MD_ZT -s workspace-write -o <scratch>\codex_last.md -
```

- 返工时续接同一会话（保留上下文）：`codex exec resume --last "<返工意见>"`，或用会话 ID 替换 `--last`。
- 只有在同一任务确实要同时改 RV_ISAR_Path 时才加 `--add-dir D:\MatlabRepo\RV_ISAR_Path`，并在任务说明中写明原因。
- **不要**使用 `--dangerously-bypass-approvals-and-sandbox`。

### 任务说明模板

```markdown
## 目标
<一句话说明要实现或修复什么>

## 背景与上下文
- 相关文件：<路径及其作用>
- 依赖的 RV_ISAR_Path / matlab_research_lib 接口：<...>

## 具体要求
1. ...

## 约束
- 只修改：<文件列表>；不要动其他文件。
- 工作区可能有用户未提交的改动，不得回退或覆盖与本任务无关的改动。
- 不要 git commit。
- 不要把 RV_ISAR_Path 中的共享函数复制到本仓库；底层算子优先复用 matlab_research_lib。

## 验收标准
- <可验证的结果，例如某个测试通过、识别结果满足某条件>
- 完成后列出修改的文件，并说明运行了哪些测试及结果。
```

## 4. 验收流程（每次 Codex 完成后必须执行）

1. **范围审查**：`git status -sb`、`git diff --stat` 与基线对比，确认没有越界修改、没有误删或覆盖用户已有改动，也没有把共享函数复制进本仓库。
2. **代码审查**：逐文件阅读 `git diff`，检查逻辑正确性、物理量与单位、矩阵维度、接口契约，以及是否违反项目规范。
3. **运行测试**（MATLAB 位于 `D:\Software\Matlab\bin\matlab.exe`；`startup_project` 会自动加载 matlab_research_lib 和 RV_ISAR_Path），优先只跑与改动相关的测试：

```powershell
matlab -batch "cd('D:\MatlabRepo\MD_ZT'); startup_project; r = runtests('tests/test_isar_ca_run_contract.m'); disp(table(r)); assert(all([r.Passed]))"
```

   - 共享 ISARpath 对接：`test_isar_micro_motion`、`test_isar_fractional_alignment`、`test_isarpath_full_chain`、`test_isarpath_native_frontend`、`test_isarpath_script_integration`。
   - 盲模式识别 / 尾点自旋：`test_isar_cone_point_kinematics`、`test_isar_tip_motion_mode_identification`、`test_isar_spin_from_tail_points`。
   - 拟合加速度 CA：`test_isar_cone_point_input_validation`、`test_isar_chain_prior_isolation`、`test_isar_golden_lo_offsets`、`test_isar_tip_to_coherent_ratio`、`test_isar_ca_run_contract`、`test_ideal_point_bench_ca`。
   - 如果 MATLAB 在沙箱内因许可证或权限原因无法启动，再改用 `BypassSandbox: true`。
4. **不合格就返工**：把具体问题（diff 位置、测试报错原文、期望行为）通过 `resume` 发回给 Codex。**最多返工 3 轮**；仍不通过就停止，向用户说明卡在哪里。

## 5. 向用户汇报

验收通过后简要说明：完成了什么、修改了哪些文件（注明属于哪个仓库）、运行了哪些测试及结果、返工了几轮、有哪些遗留风险或需要用户决定的问题。
