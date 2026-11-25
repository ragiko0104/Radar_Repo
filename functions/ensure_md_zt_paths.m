function paths = ensure_md_zt_paths(verbose)
%ENSURE_MD_ZT_PATHS Function-safe path bootstrap for MD_ZT.
%   Adds MD_ZT folders, matlab_research_lib/lib, and RV_ISAR_Path (code/
%   and functions/) without polluting caller static function workspaces.

if nargin < 1 || isempty(verbose)
    verbose = false;
end

projectRoot = fileparts(fileparts(mfilename('fullpath')));

% 1. Add MD_ZT folders
subdirs = {'code', 'functions', 'experiments', 'tests'};
for idx = 1:numel(subdirs)
    sub = fullfile(projectRoot, subdirs{idx});
    if exist(sub, 'dir') == 7
        addpath(genpath(sub));
    end
end

% Ensure local output folders exist
for folderName = {'results', 'cache'}
    folderPath = fullfile(projectRoot, folderName{1});
    if exist(folderPath, 'dir') ~= 7
        mkdir(folderPath);
    end
end

% Load optional local_config.m into isolated function workspace
localCfg = struct();
localConfigFile = fullfile(projectRoot, 'local_config.m');
if exist(localConfigFile, 'file') == 2
    localCfg = load_local_config(localConfigFile);
end

% 2. Resolve matlab_research_lib/lib
libCandidates = {};
if isfield(localCfg, 'MATLAB_RESEARCH_LIB_ROOT') && ~isempty(localCfg.MATLAB_RESEARCH_LIB_ROOT)
    libCandidates{end+1} = char(string(localCfg.MATLAB_RESEARCH_LIB_ROOT));
end
envLib = getenv('MATLAB_RESEARCH_LIB_ROOT');
if ~isempty(envLib)
    libCandidates{end+1} = envLib;
end
libCandidates{end+1} = fullfile(projectRoot, 'external', 'matlab_research_lib');
libCandidates{end+1} = fullfile(fileparts(projectRoot), 'matlab_research_lib');
libCandidates{end+1} = 'D:\MatlabProj\matlab_research_lib';

libRoot = '';
libDir = '';
for idx = 1:numel(libCandidates)
    candLib = fullfile(libCandidates{idx}, 'lib');
    if exist(candLib, 'dir') == 7
        libRoot = libCandidates{idx};
        libDir = candLib;
        break;
    end
end
if isempty(libDir)
    error('MD_ZT:MissingResearchLib', ...
        ['matlab_research_lib/lib not found. Set MATLAB_RESEARCH_LIB_ROOT ' ...
        'or define it in local_config.m (copy from local_config_template.m).']);
end
addpath(genpath(libDir));

% 3. Resolve shared RV_ISAR_Path repository
rvCandidates = {};
if isfield(localCfg, 'RV_ISAR_PATH_ROOT') && ~isempty(localCfg.RV_ISAR_PATH_ROOT)
    rvCandidates{end+1} = char(string(localCfg.RV_ISAR_PATH_ROOT));
end
envRv = getenv('RV_ISAR_PATH_ROOT');
if ~isempty(envRv)
    rvCandidates{end+1} = envRv;
end
rvCandidates{end+1} = fullfile(fileparts(projectRoot), 'RV_ISAR_Path');
rvCandidates{end+1} = 'D:\MatlabRepo\RV_ISAR_Path';

rvRoot = '';
for idx = 1:numel(rvCandidates)
    if exist(fullfile(rvCandidates{idx}, 'code', 'multiuav_isarpath.m'), 'file') == 2
        rvRoot = rvCandidates{idx};
        break;
    end
end
if isempty(rvRoot)
    error('MD_ZT:MissingRvIsarPath', ...
        ['RV_ISAR_Path repository not found. Place RV_ISAR_Path beside MD_ZT, ' ...
        'set RV_ISAR_PATH_ROOT, or define it in local_config.m.']);
end
addpath(genpath(fullfile(rvRoot, 'code')));
addpath(genpath(fullfile(rvRoot, 'functions')));
if exist(fullfile(rvRoot, 'config'), 'dir') == 7
    addpath(genpath(fullfile(rvRoot, 'config')));
end

paths = struct( ...
    'projectRoot', projectRoot, ...
    'rvIsarPathRoot', rvRoot, ...
    'matlabResearchLibRoot', libRoot, ...
    'matlabResearchLibDir', libDir);

if verbose
    fprintf('MD_ZT startup complete.\n');
    fprintf('Project root: %s\n', projectRoot);
    fprintf('RV_ISAR_Path root: %s\n', rvRoot);
    fprintf('MATLAB research library: %s\n', libDir);
end
end

function cfg = load_local_config(localConfigFile)
MATLAB_RESEARCH_LIB_ROOT = ''; %#ok<NASGU>
RV_ISAR_PATH_ROOT = ''; %#ok<NASGU>
DATA_ROOT = ''; %#ok<NASGU>
RESULT_ROOT = ''; %#ok<NASGU>
run(localConfigFile);
cfg = struct( ...
    'MATLAB_RESEARCH_LIB_ROOT', MATLAB_RESEARCH_LIB_ROOT, ...
    'RV_ISAR_PATH_ROOT', RV_ISAR_PATH_ROOT, ...
    'DATA_ROOT', DATA_ROOT, ...
    'RESULT_ROOT', RESULT_ROOT);
end
