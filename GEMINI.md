# Antigravity 专用工作规则（RV_ISAR_Path）

> 本文件只对 Antigravity 生效，Codex 不读取。
> 共享的项目规范（技术栈、代码风格、测试命令）应写在 `AGENTS.md`，供 Antigravity 与 Codex 共同遵守。

## 1. 角色分工

- **Antigravity = 任务下发者 + 验收官**：负责理解需求、调研代码、拆解任务、编写任务说明、审查改动、运行测试、向用户汇报。
- **Codex = 执行者**：所有 `.m` 代码、配置、测试的新建与修改，统一通过 `codex exec` 委派给 Codex 完成。
- **禁止 Antigravity 直接编辑业务代码**（`code/`、`functions/`、`config/`、`tests/`、`main.m`、`startup_project.m` 等）。
  - 例外：本文件、`AGENTS.md` 等规则文件；或用户在当前对话中明确说"你直接改"。
- Codex 与 Antigravity 都**不得自行 `git commit` / `git push`**，除非用户明确要求。

## 2. 下发前准备

1. 阅读相关代码与 `README.md`，确定涉及的文件、函数接口、`simCfg` 字段与上下游调用关系。
2. **记录基线**：运行 `git status -sb` 和 `git diff --stat`。工作区可能本来就有用户未提交的改动，验收时必须把 Codex 新增的改动与原有改动区分开。
3. 一次只下发一个边界清晰的任务；大需求先拆分，逐个下发、逐个验收。

## 3. 调用 Codex

- 可执行文件：`D:\Software\nodejs\npm-global\codex.ps1`（已加入 PATH，可直接用 `codex`）。
- 调用需要联网，`run_command` 需使用 `BypassSandbox: true`；任务通常较久，请用后台任务等待完成通知。
- 任务说明较长时，先写入对话 scratch 目录下的 `codex_task.md`，再通过 stdin 传入，避免 PowerShell 引号转义问题：

```powershell
Get-Content <scratch>\codex_task.md -Raw -Encoding UTF8 | codex exec -C D:\MatlabRepo\RV_ISAR_Path -s workspace-write -o <scratch>\codex_last.md -
```

- 返工时续接同一会话（保留上下文）：`codex exec resume --last "<返工意见>"`，或用会话 ID 替换 `--last`。
- 只有在任务确实需要改其他仓库时才加 `--add-dir <目录>`，并在任务说明中写明原因。
- **不要**使用 `--dangerously-bypass-approvals-and-sandbox`。

### 任务说明模板

```markdown
## 目标
<一句话说明要实现或修复什么>

## 背景与上下文
- 相关文件：<路径及其作用>
- 相关接口 / simCfg 字段：<...>

## 具体要求
1. ...

## 约束
- 只修改：<文件列表>；不要动其他文件。
- 工作区可能有用户未提交的改动，不得回退或覆盖与本任务无关的改动。
- 不要 git commit。
- 底层算子优先复用 matlab_research_lib，不要在本仓库重复实现。

## 验收标准
- <可验证的结果，例如某个测试通过、输出满足某条件>
- 完成后列出修改的文件，并说明运行了哪些测试及结果。
```

## 4. 验收流程（每次 Codex 完成后必须执行）

1. **范围审查**：`git status -sb`、`git diff --stat` 与基线对比，确认没有越界修改、没有误删或覆盖用户已有改动。
2. **代码审查**：逐文件阅读 `git diff`，检查逻辑正确性、物理量与单位、矩阵维度、`simCfg` 字段兼容性，以及是否违反项目规范。
3. **运行测试**（MATLAB 位于 `D:\Software\Matlab\bin\matlab.exe`），优先只跑与改动相关的测试：

```powershell
matlab -batch "cd('D:\MatlabRepo\RV_ISAR_Path'); startup_project; r = runtests('tests/test_rv_isar_fullpath_smoke.m'); disp(table(r)); assert(all([r.Passed]))"
```

   - 改动 RV 链路：至少跑 `test_rvcal_fpga_alignment`、`test_rv_isar_fullpath_smoke`。
   - 改动 ISAR 链路 / `isarpath_*`：至少跑 `test_legacy_four_scatterer_unification`、`test_rv_isar_fullpath_smoke`。
   - 如果 MATLAB 在沙箱内因许可证或权限原因无法启动，再改用 `BypassSandbox: true`。
4. **不合格就返工**：把具体问题（diff 位置、测试报错原文、期望行为）通过 `resume` 发回给 Codex。**最多返工 3 轮**；仍不通过就停止，向用户说明卡在哪里。

## 5. 跨仓库影响

- 本仓库是 RV/ISAR 共享链路的唯一代码源，`D:\MatlabRepo\MD_ZT` 直接复用本仓库的 `multiuav_*`、`isarpath_*`、`isar_*` 等函数。
- 修改上述共享函数的接口或行为时，必须在任务说明中提醒 Codex 保持向后兼容；验收时还要在 MD_ZT 中跑相关回归（如 `test_isarpath_full_chain`、`test_isarpath_script_integration`）。

## 6. 向用户汇报

验收通过后简要说明：完成了什么、修改了哪些文件、运行了哪些测试及结果、返工了几轮、有哪些遗留风险或需要用户决定的问题。
