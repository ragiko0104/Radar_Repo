function savedPath = save_output_artifact(artifact, artifactKind, saveOpts, defaultBaseName)
%SAVE_OUTPUT_ARTIFACT 统一保存图形对象或评估表格/结构体到指定目录。
%
% 功能说明：
%   为所有 plot_* 与 assess_* 函数提供一致的可选落盘行为：
%   - 当 saveOpts.enable == false 或 saveOpts 为空时直接跳过并返回 ""；
%   - 当 saveOpts.enable == true 时自动创建目标目录（如不存在）并按指定格式保存。
%
% 输入参数：
%   artifact        - 待保存对象：
%                       artifactKind="figure" 时为 1×1 Figure 句柄；
%                       artifactKind="table"  时为 MATLAB table。
%   artifactKind    - 1×1 string/char，"figure" 或 "table"
%   saveOpts        - (可选) 1×1 struct，保存配置结构体：
%     .enable         1×1 logical，是否保存（默认 false）
%     .dir            char/string，保存目录路径（enable==true 时必填）
%     .fileName       (可选) char/string，自定义文件名（可带或不带扩展名）
%     .format         (可选) char/string：
%                       figure 默认 'png'（支持 'png'/'pdf'/'jpg'/'eps'/'fig'）；
%                       table  默认 'csv'（支持 'csv'/'xlsx'/'mat'）
%     .dpi            (可选) 1×1 double，位图导出分辨率 DPI（默认 300）
%   defaultBaseName - 1×1 string/char，未指定 fileName 时使用的默认主文件名
%
% 输出参数：
%   savedPath       - 1×1 string，实际写入的文件完整路径（未保存时返回 ""）

savedPath = "";
if nargin < 3 || isempty(saveOpts) || ~isstruct(saveOpts)
    return;
end
if ~isfield(saveOpts, 'enable') || ~logical(saveOpts.enable)
    return;
end

if ~isfield(saveOpts, 'dir') || isempty(saveOpts.dir)
    error('RV_ISAR_Path:MissingSaveDir', ...
        'saveOpts.enable=true 时必须提供非空的 saveOpts.dir 路径。');
end
outDir = char(string(saveOpts.dir));
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

kind = lower(string(artifactKind));
if kind == "figure"
    defaultFmt = 'png';
elseif kind == "table"
    defaultFmt = 'csv';
else
    error('RV_ISAR_Path:UnsupportedArtifactKind', ...
        '不支持的 artifactKind="%s"，仅支持 "figure" 或 "table"。', kind);
end

fmt = defaultFmt;
if isfield(saveOpts, 'format') && ~isempty(saveOpts.format)
    fmt = lower(char(string(saveOpts.format)));
    fmt = regexprep(fmt, '^\.', '');
end

if isfield(saveOpts, 'fileName') && ~isempty(saveOpts.fileName)
    rawName = char(string(saveOpts.fileName));
    [~, baseName, ext] = fileparts(rawName);
    if isempty(ext)
        fullPath = fullfile(outDir, [baseName, '.', fmt]);
    else
        fullPath = fullfile(outDir, rawName);
        fmt = lower(regexprep(ext, '^\.', ''));
    end
else
    fullPath = fullfile(outDir, [char(string(defaultBaseName)), '.', fmt]);
end

if kind == "figure"
    dpi = 300;
    if isfield(saveOpts, 'dpi') && ~isempty(saveOpts.dpi)
        dpi = double(saveOpts.dpi);
    end
    if strcmp(fmt, 'fig')
        savefig(artifact, fullPath);
    elseif exist('exportgraphics', 'file') == 2 && ismember(fmt, {'png', 'jpg', 'jpeg', 'tif', 'tiff', 'pdf', 'eps'})
        exportgraphics(artifact, fullPath, 'Resolution', dpi);
    else
        saveas(artifact, fullPath);
    end
else
    if strcmp(fmt, 'mat')
        resultTable = artifact; %#ok<NASGU>
        save(fullPath, 'resultTable');
    else
        writetable(artifact, fullPath);
    end
end

savedPath = string(fullPath);
end
