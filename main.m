%% main.m
% RV_ISAR_Path 外部参数配置与纯计算函数化运行总入口。
%
% 使用方式：
%   1. 在下方修改 simCfg 的运行模式（"rv" / "isar" / "both"）、目标参数、
%      散射体模型、波形/本振偏置/RCMC 后处理参数，然后运行 main；
%   2. 运行完成后，工作区保留 result 与 scene：
%        - 需要评估指标与导出表格：运行 code/assess_rv_isar_results.m
%        - 需要绘制中间量/成像图/动图与导出图片：运行 code/plot_rv_isar_results.m
%   3. 若需运行紧凑 UAV 场景 A–E，可将 targetModel = "uavScatterer" 后运行 main。

startup_project;

if ~exist('targetModel', 'var') || isempty(targetModel)
    targetModel = "fullpath";
end

switch lower(string(targetModel))
    case "fullpath"
        %% 1. 加载独立默认配置（如工作区已有 simCfg 则直接在其上运行）
        if ~exist('simCfg', 'var') || isempty(simCfg)
            simCfg = rv_isar_default_config();

            %% 2. 链路选择与核心计算开关（纯计算，不自动弹窗出图）
            simCfg.switches.runMode           = "both";      % "rv" | "isar" | "both"
            simCfg.switches.rvLfmSpectrumMode = 'floating';  % 'floating' | 'fixedPoint' | 'idealTime'

            %% 3. 目标个数、运动参数与散射体模型设置（统一调用 lib 预设函数）
            % simCfg.targets.targetIndices     = [1 2];
            % simCfg.targets.motionLosAngleDeg = 45;
            simCfg.scatterer.modelTypes      = "legacy_ZXT4"; % legacy_ZXT4 / legacy_plane4 / rigid_cross4 / ideal_point / quad_small ...
            % simCfg.scatterer.overrides       = struct();                               % 可选：局部覆盖预设参数，如 struct('radiusM', 0.5)

            %% 4. RV / ISAR 波形、本振拍频偏置、回波性能选项与 RCMC/后处理参数按需覆盖
            % simCfg.rv.Na                                    = 1;
            % simCfg.isar.Na                                  = 32;
            % simCfg.isar.echoMode                            = "stopAndGo";        % "fastTimeExact" | "stopAndGo"
            % simCfg.isar.frequencyBiasHz                     = [20e6, 45e6];       % 正=本振在目标左侧，负=右侧
            simCfg.isar.postprocess.rcmcSearchRadiusBins    = 80;                 % RCMC 最大搜索位移半径 (±bins)
            % simCfg.isar.postprocess.phaseCompensationMode   = "VARIANCE_UNIT";    % "VARIANCE_UNIT" | "MODEL_CV" | "MODEL_CA" | "NONE"
            %
            % 旋转与微动模式设置（按需取消注释）：
            % 1) 关闭微动，仅绕自定义体坐标系转轴匀速转动（如绕 +z 轴 1 Hz）：
            simCfg.isar.bodyRotationMode          = "constantBodyRate";
            simCfg.isar.omegaBody_B               = [0, 0, pi/2];       % [wx, wy, wz] (rad/s)
            % 2) 开启锥旋（coneSpin）/ 摆动（swingSpin）/ 复合锥摆旋（coneSwingSpin）：
            % simCfg.isar.bodyRotationMode          = "coneSwingSpin";    % "spin" | "coneSpin" | "swingSpin" | "coneSwingSpin"
            % simCfg.isar.spinFrequencyHz           = 1;                  % 自旋频率 (Hz)
            % simCfg.isar.precessionFrequencyHz     = 1;                  % 锥旋频率 (Hz)
            % simCfg.isar.precessionHalfAngleDeg    = 10;                 % 锥旋半顶角 (deg)
            % simCfg.isar.swingFrequencyHz          = 1;                  % 摆动频率 (Hz)
            % simCfg.isar.swingAngleMinDeg          = -5;                 % 摆动最小角 (deg)
            % simCfg.isar.swingAngleMaxDeg          = 5;                  % 摆动最大角 (deg)
            %
            % 回波分组与分块选项（默认均关闭以获得最高速度；一般在仿真 UAV 等数百点散射体时再启用分块控内存）：
            % simCfg.rv.returnGroups                          = false;              % true 时保留 body/rotor 等分组回波
            % simCfg.isar.returnGroups                        = false;              % true 时保留各散射体分组回波
            % simCfg.rv.fastTimeChunkSize                     = [];                 % 如仿 UAV 时可设 4096
            % simCfg.rv.scattererChunkSize                    = [];                 % 如仿 UAV 时可设 96
            % simCfg.isar.fastTimeChunkSize                   = [];                 % 如仿 UAV 时可设 4096
            % simCfg.isar.scattererChunkSize                  = [];                 % 如仿 UAV 时可设 96
        end

        %% 5. 调用纯计算全链路总入口
        [result, scene] = rv_isar_run(simCfg);

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
