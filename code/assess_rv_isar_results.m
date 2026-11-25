%% assess_rv_isar_results.m
% RV_ISAR_Path 独立评估总入口脚本（计算与评估彻底解耦）。
%
% 使用前提：
%   工作区已有通过 main.m 或 [result, scene] = rv_isar_run(simCfg) 产出的 result。
% 使用方式：
%   1. 在顶部配置结果保存路径 saveDir、是否保存 enableSave、目标范围与评估项开关；
%   2. 可整脚本运行，也可按分节（Ctrl+Enter）单独执行某一项评估。

if isempty(which('assess_rv_measurement'))
    run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup_project.m'));
end

%% ==================== 顶部统一保存与评估配置入口 ====================
projectRoot = fileparts(fileparts(mfilename('fullpath')));
date = datetime("today");
id = 1;
saveDir     = fullfile(projectRoot, 'results', 'assess'); % 表格保存根目录（可自行修改）
enableSave  = false;                                      % true: 自动保存表格；false: 仅计算显示
saveFormat  = 'csv';                                      % 'csv' | 'xlsx' | 'mat'

assessCfg = struct( ...
    'targetIndices',     [], ...        % 空表示全部目标，可指定如 [1 2]
    'rvPulseIndex',      1, ...         % RV 评估脉冲索引
    'expectedPeakCount', 4, ...         % ISAR 散射点评估期望亮点数
    'resolutionMethod',  "linear", ...  % 3dB 宽度计算方法
    'printSummary',      true);         % 是否向命令行打印各评估表

saveOptsBase = struct('enable', enableSave, 'dir', saveDir, 'format', saveFormat);

% 各项独立评估开关（按需设为 true / false）
doAssessRv             = true;
doAssessIsarDispersion = true;
doAssessIsarRcmc       = true;
doAssessIsarImage      = true;
doAssessIsarScatterers = true;
doSaveCombinedIsarTbl  = true;

%% ==================== 1. RV 测距测速精度评估 ====================
if doAssessRv && exist('result', 'var') && isfield(result, 'ranRV') && result.ranRV
    rvCfg = struct( ...
        'pulseIndex',    assessCfg.rvPulseIndex, ...
        'targetIndices', assessCfg.targetIndices, ...
        'printSummary',  assessCfg.printSummary);
    [rvAssessDetail, rvAssessTable, rvSavedPath] = assess_rv_measurement( ...
        result.rv, rvCfg, saveOptsBase);
end

%% ==================== 2. ISAR 色散补偿前后 HRRP 质量评估 ====================
if doAssessIsarDispersion && exist('result', 'var') && isfield(result, 'ranISAR') && result.ranISAR
    [dispDetail, dispTable, dispSavedPath] = assess_isar_dispersion( ...
        result.isar, assessCfg, saveOptsBase);
end

%% ==================== 3. ISAR RCMC 包络对齐质量评估 ====================
if doAssessIsarRcmc && exist('result', 'var') && isfield(result, 'ranISAR') && result.ranISAR
    [rcmcDetail, rcmcTable, rcmcSavedPath] = assess_isar_rcmc( ...
        result.isar, assessCfg, saveOptsBase);
end

%% ==================== 4. ISAR 图像聚焦质量评估（熵/PSLR/3dB 能量比）====================
if doAssessIsarImage && exist('result', 'var') && isfield(result, 'ranISAR') && result.ranISAR
    [imgAssessDetail, imgAssessTable, imgSavedPath] = assess_isar_image( ...
        result.isar, assessCfg, saveOptsBase);
end

%% ==================== 5. ISAR 散射点二维 3dB 分辨率评估 ====================
if doAssessIsarScatterers && exist('result', 'var') && isfield(result, 'ranISAR') && result.ranISAR
    [scatDetail, scatTable, scatSavedPath] = assess_isar_scatterers( ...
        result.isar, assessCfg, saveOptsBase);
end

%% ==================== 6. 汇总全部 ISAR 评估指标为总表 ====================
if exist('dispTable', 'var') && exist('rcmcTable', 'var') && ...
        exist('imgAssessTable', 'var') && exist('scatTable', 'var')
    isarAssessSummary = join(join(join(dispTable, rcmcTable, 'Keys', 'Target'), ...
        imgAssessTable, 'Keys', 'Target'), scatTable, 'Keys', 'Target');
    if assessCfg.printSummary
        fprintf('\n=== ISAR 全项质量评估总表 ===\n');
        disp(isarAssessSummary);
    end
    if doSaveCombinedIsarTbl
        isarSummarySavedPath = save_output_artifact( ...
            isarAssessSummary, "table", saveOptsBase, 'assess_isar_summary');
    end
end
