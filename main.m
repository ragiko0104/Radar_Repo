%% main.m
% Top-level entry for RV_ISAR_Path (lib-backed RV + ISAR workflows).
%
% Optional caller workspace overrides before running main:
%   targetModel   - "fullpath" (default, runs code/multiuav_iasr_fullpath.m)
%                   or "uavScatterer" (runs compact uavmodel A-E scenarios)
%   uavScenario   - "A", "B", "C", "D", "E", or "ALL" (uavScatterer mode)
%   uavEchoMode   - "stopAndGo" or "fastTimeExact" (uavScatterer mode)
%   run_RVpath    - logical override for fullpath RV branch
%   run_ISARpath  - logical override for fullpath ISAR branch
%   targetModelTypes - "uav", "legacyDisk", "rigidPointCloud", "point",
%                      or "projectedRange" (fullpath mode)

startup_project;

if ~exist('targetModel', 'var') || isempty(targetModel)
    targetModel = "fullpath";
end

switch lower(string(targetModel))
    case "fullpath"
        run(fullfile(project_root, 'code', 'multiuav_iasr_fullpath.m'));

    case "uavscatterer"
        if isempty(which('uavmodel.createModel'))
            error('RV_ISAR_Path:MissingUavModel', ...
                ['The UAV model library is not available. Configure ' ...
                'MATLAB_RESEARCH_LIB_ROOT or local_config.m.']);
        end
        if ~exist('uavScenario', 'var') || isempty(uavScenario)
            uavScenario = "A";
        end
        if ~exist('uavEchoMode', 'var') || isempty(uavEchoMode)
            uavEchoMode = "stopAndGo";
        end
        uav_result = run_uav_existing_pipeline(uavScenario, ...
            struct('EchoMode', uavEchoMode, 'ShowFigures', true));

    otherwise
        error('RV_ISAR_Path:UnsupportedTargetModel', ...
            'Unsupported targetModel "%s". Use "fullpath" or "uavScatterer".', ...
            string(targetModel));
end
