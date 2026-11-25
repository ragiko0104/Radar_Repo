%% startup_project.m
% Project startup script for RV_ISAR_Path (Radar_Repo).
% Adds project folders and resolves matlab_research_lib/lib cleanly.

% Get project root directory
project_root = fileparts(mfilename('fullpath'));

% Add existing project folders to MATLAB path
project_subdirs = {'config', 'code', 'functions', 'tests'};
for idx = 1:numel(project_subdirs)
    subdir_path = fullfile(project_root, project_subdirs{idx});
    if exist(subdir_path, 'dir') == 7
        addpath(genpath(subdir_path));
    end
end

% Create local result folders if they do not exist
result_dir = fullfile(project_root, 'results');
if exist(result_dir, 'dir') ~= 7
    mkdir(result_dir);
end

cache_dir = fullfile(project_root, 'cache');
if exist(cache_dir, 'dir') ~= 7
    mkdir(cache_dir);
end

% Load optional local machine-specific configuration
local_config_file = fullfile(project_root, 'local_config.m');
if exist(local_config_file, 'file') == 2
    run(local_config_file);
end

% Resolve matlab_research_lib root in priority order:
%   1. local_config.m / caller workspace variable MATLAB_RESEARCH_LIB_ROOT
%   2. Environment variable MATLAB_RESEARCH_LIB_ROOT
%   3. Project-local external/matlab_research_lib
%   4. Sibling folder ../matlab_research_lib
%   5. Local workstation default D:\MatlabProj\matlab_research_lib
lib_candidates = {};
if exist('MATLAB_RESEARCH_LIB_ROOT', 'var') && ~isempty(MATLAB_RESEARCH_LIB_ROOT)
    lib_candidates{end+1} = char(string(MATLAB_RESEARCH_LIB_ROOT));
end
env_lib_root = getenv('MATLAB_RESEARCH_LIB_ROOT');
if ~isempty(env_lib_root)
    lib_candidates{end+1} = env_lib_root;
end
lib_candidates{end+1} = fullfile(project_root, 'external', 'matlab_research_lib');
lib_candidates{end+1} = fullfile(fileparts(project_root), 'matlab_research_lib');
lib_candidates{end+1} = 'D:\MatlabProj\matlab_research_lib';

matlab_research_lib_root = '';
matlab_research_lib_dir = '';
for idx = 1:numel(lib_candidates)
    candidate_root = lib_candidates{idx};
    candidate_lib = fullfile(candidate_root, 'lib');
    if exist(candidate_lib, 'dir') == 7
        matlab_research_lib_root = candidate_root;
        matlab_research_lib_dir = candidate_lib;
        break;
    end
end

if ~isempty(matlab_research_lib_dir)
    addpath(genpath(matlab_research_lib_dir));
else
    error('RV_ISAR_Path:MissingResearchLib', ...
        ['matlab_research_lib/lib not found. Set environment variable ' ...
        'MATLAB_RESEARCH_LIB_ROOT or define MATLAB_RESEARCH_LIB_ROOT in ' ...
        'local_config.m (copy from local_config_template.m).']);
end

% Display project information
fprintf('RV_ISAR_Path startup complete.\n');
fprintf('Project root: %s\n', project_root);
fprintf('MATLAB research library: %s\n', matlab_research_lib_dir);
