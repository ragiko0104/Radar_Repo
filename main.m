%% main.m
% RV_ISAR_Path 外部参数配置与运行总入口（独立于 fullpath 核心链路脚本）。
%
% 使用方式：
%   1. 在下方直接修改 simCfg 的运行模式（"rv" / "isar" / "both"）、目标参数、
%      散射体模型、波形参数与绘图开关，然后直接运行 main；
%   2. 或在工作区预先构造 simCfg 后调用 run("code/multiuav_iasr_fullpath.m")；
%   3. 若需运行紧凑 UAV 场景 A–E，可将 targetModel = "uavScatterer" 后运行 main。

startup_project;

if ~exist('targetModel', 'var') || isempty(targetModel)
    targetModel = "fullpath";
end

switch lower(string(targetModel))
    case "fullpath"
        %% 1. 加载独立默认配置（如工作区已有 simCfg 则在其基础上补齐默认值）
        if ~exist('simCfg', 'var') || isempty(simCfg)
            simCfg = rv_isar_default_config();

            %% 2. 链路选择与开关（外部控制只跑 RV_path / 只跑 ISAR_path / 都跑）
            simCfg.switches.runMode = "both";          % "rv" | "isar" | "both"
            simCfg.switches.uavPlotGeometryEnabled  = true;
            simCfg.switches.rvShowEchoTimeFrequency = true;
            simCfg.switches.rvCfarPlotEnabled       = true;
            simCfg.switches.isarPlotEnabled         = true;
            simCfg.switches.isarQualityEnabled      = true;
            simCfg.switches.isarQualityPrintEnabled = true;

            %% 3. 目标个数、运动参数与散射体模型设置
            % 提示：可直接修改 R0/v_Mach/a_g，或通过 targetIndices 选取子集（如 [1 2]）
            % simCfg.targets.targetIndices = [1 2];
            % simCfg.targets.motionLosAngleDeg = 45;
            % simCfg.scatterer.modelTypes = "legacyDisk"; % "uav" | "legacyDisk" | "rigidPointCloud" | "point" | "projectedRange"

            %% 4. RV / ISAR 波形参数按需覆盖
            % simCfg.rv.Na = 1;
            % simCfg.isar.Na = 32;
            % simCfg.isar.echoMode = "fastTimeExact";  % ISAR 可选 "stopAndGo"（快时间冻结）
        end

        %% 5. 运行全链路执行脚本
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
