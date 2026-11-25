%% local_config_template.m
% Copy this file to local_config.m and modify it for each computer.
MATLAB_RESEARCH_LIB_ROOT = 'D:\MatlabProj\matlab_research_lib';
RV_ISAR_PATH_ROOT = 'D:\MatlabRepo\RV_ISAR_Path';
DATA_ROOT = 'D:\RadarData';
RESULT_ROOT = fullfile(fileparts(mfilename('fullpath')), 'results');

fprintf('MATLAB_RESEARCH_LIB_ROOT: %s\n', MATLAB_RESEARCH_LIB_ROOT);
fprintf('RV_ISAR_PATH_ROOT: %s\n', RV_ISAR_PATH_ROOT);
fprintf('DATA_ROOT: %s\n', DATA_ROOT);
fprintf('RESULT_ROOT: %s\n', RESULT_ROOT);
