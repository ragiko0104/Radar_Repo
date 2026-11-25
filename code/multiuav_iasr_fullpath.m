%MULTIUAV_IASR_FULLPATH RV + ISAR 全链路纯计算兼容薄脚本（调用 rv_isar_build_scene 与 rv_isar_run）。
% 推荐新代码通过 main.m 或直接调用 [result, scene] = rv_isar_run(simCfg)。
% 绘图请调用 code/plot_rv_isar_results.m，评估请调用 code/assess_rv_isar_results.m。

if isempty(which('radar.echo_generate_target_lfm')) || isempty(which('rv_isar_default_config'))
    run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup_project.m'));
end

%% 加载并合并外部配置与兼容工作区覆盖项
if exist('simCfg', 'var') && ~isempty(simCfg)
    simCfg = rv_isar_default_config(simCfg);
else
    simCfg = rv_isar_default_config();
end

c = simCfg.constants.c;
if exist('run_RVpath', 'var') && ~isempty(run_RVpath) && ...
        exist('run_ISARpath', 'var') && ~isempty(run_ISARpath)
    if logical(run_RVpath) && logical(run_ISARpath)
        simCfg.switches.runMode = "both";
    elseif logical(run_RVpath)
        simCfg.switches.runMode = "rv";
    elseif logical(run_ISARpath)
        simCfg.switches.runMode = "isar";
    end
elseif exist('run_RVpath', 'var') && ~isempty(run_RVpath)
    if logical(run_RVpath), simCfg.switches.runMode = "both"; else, simCfg.switches.runMode = "isar"; end
elseif exist('run_ISARpath', 'var') && ~isempty(run_ISARpath)
    if logical(run_ISARpath), simCfg.switches.runMode = "both"; else, simCfg.switches.runMode = "rv"; end
end

if exist('rvLfmSpectrumMode', 'var') && ~isempty(rvLfmSpectrumMode)
    simCfg.switches.rvLfmSpectrumMode = rvLfmSpectrumMode;
end

if exist('FULLPATH_OVERRIDES', 'var') && isstruct(FULLPATH_OVERRIDES)
    if isfield(FULLPATH_OVERRIDES, 'isarNa') && ~isempty(FULLPATH_OVERRIDES.isarNa)
        simCfg.isar.Na = double(FULLPATH_OVERRIDES.isarNa);
    end
    if isfield(FULLPATH_OVERRIDES, 'isarEchoMode') && ~isempty(FULLPATH_OVERRIDES.isarEchoMode)
        simCfg.isar.echoMode = string(FULLPATH_OVERRIDES.isarEchoMode);
    end
    if isfield(FULLPATH_OVERRIDES, 'targetIndices') && ~isempty(FULLPATH_OVERRIDES.targetIndices)
        simCfg.targets.targetIndices = FULLPATH_OVERRIDES.targetIndices;
    end
end
if exist('targetModelTypes', 'var') && ~isempty(targetModelTypes)
    simCfg.scatterer.modelTypes = targetModelTypes;
end

if exist('ISAR_MOTION_OPTIONS', 'var') && ~isempty(ISAR_MOTION_OPTIONS)
    motionOptionNames__ = {'Mode', 'SpinFrequencyHz', ...
        'PrecessionFrequencyHz', 'PrecessionHalfAngleDeg'};
    motionConfigNames__ = {'bodyRotationMode', 'spinFrequencyHz', ...
        'precessionFrequencyHz', 'precessionHalfAngleDeg'};
    if ~(isstruct(ISAR_MOTION_OPTIONS) && isscalar(ISAR_MOTION_OPTIONS)) || ...
            any(~ismember(fieldnames(ISAR_MOTION_OPTIONS), motionOptionNames__))
        error('RadarProject:ISARMotion:InvalidOptions', ...
            'ISAR_MOTION_OPTIONS contains an unsupported field.');
    end
    for k__ = 1:numel(motionOptionNames__)
        if isfield(ISAR_MOTION_OPTIONS, motionOptionNames__{k__})
            simCfg.isar.(motionConfigNames__{k__}) = ISAR_MOTION_OPTIONS.(motionOptionNames__{k__});
        end
    end
    clear motionOptionNames__ motionConfigNames__ k__
end
if exist('ISAR_POST_OPTIONS', 'var') && ~isempty(ISAR_POST_OPTIONS)
    postNames__ = fieldnames(ISAR_POST_OPTIONS);
    for k__ = 1:numel(postNames__)
        simCfg.isar.postprocess.(postNames__{k__}) = ISAR_POST_OPTIONS.(postNames__{k__});
    end
    clear postNames__ k__
end
if exist('ISAR_PRIOR_ACCELERATION_MPS2', 'var') && ~isempty(ISAR_PRIOR_ACCELERATION_MPS2)
    simCfg.isar.priorAccelerationMps2 = ISAR_PRIOR_ACCELERATION_MPS2;
end

cfg_RVpath = simCfg.rv;
cfg_ISARpath = simCfg.isar;
cfg_ISARpath.propagationSpeed = c;
cfg_ISARpath.B      = abs(cfg_ISARpath.fH - cfg_ISARpath.fL);
cfg_ISARpath.kr     = (cfg_ISARpath.fH - cfg_ISARpath.fL) / cfg_ISARpath.Tr;
cfg_ISARpath.fc     = (cfg_ISARpath.fH + cfg_ISARpath.fL) / 2;
cfg_ISARpath.lambda = cfg_ISARpath.propagationSpeed / cfg_ISARpath.fc;
cfg_ISARpath.omegaBody_B = [0, 2*pi*cfg_ISARpath.spinFrequencyHz, 0];

%% 其余仿真参数定义
scene = rv_isar_build_scene(simCfg);
if exist('run_RVpath', 'var') && ~isempty(run_RVpath)
    scene.switches.runRV = logical(run_RVpath);
end
if exist('run_ISARpath', 'var') && ~isempty(run_ISARpath)
    scene.switches.runISAR = logical(run_ISARpath);
end
if ~scene.switches.runRV && ~scene.switches.runISAR
    warning('RV_ISAR_Path:NoPathSelected', ...
        'run_RVpath 与 run_ISARpath 均为 false，本次未执行任何链路。');
end

[rvIsarResult, scene] = rv_isar_run(scene);

simCfg                      = scene.cfg;
cfg_RVpath                  = scene.cfgRV;
cfg_ISARpath                = scene.cfgISAR;
radarState                  = scene.radarState;
scatterer                   = scene.scatterer;
targets                     = scene.targets;
n_Targets                   = scene.nTargets;
R0                          = scene.truth.R0;
v_Mach                      = scene.truth.v_Mach;
a_g                         = scene.truth.a_g;
v_radial0                   = scene.truth.v_radial0;
a_radial                    = scene.truth.a_radial;
initialLosMotionDiagnostics = scene.truth.motionDiagnostics;
R_isar                      = rvIsarResult.priorMotion.rangeM;
v_isar                      = rvIsarResult.priorMotion.velocityMps;

if rvIsarResult.ranRV
    rvEcho             = rvIsarResult.rv.echo;
    rvPulseCompression = rvIsarResult.rv.pulseCompression;
    rvCfar             = rvIsarResult.rv.cfar;
    rvMeasurement      = rvIsarResult.rv.measurement;
end
if rvIsarResult.ranISAR
    targetsISAR            = rvIsarResult.isar.targetsISAR;
    isarMotion             = rvIsarResult.isar.motion;
    isarCrossRangeGeometry = rvIsarResult.isar.crossRangeGeometry;
    isarScattererState     = rvIsarResult.isar.scattererState;
    isarEcho               = rvIsarResult.isar.echo;
    isarEchoInfo           = rvIsarResult.isar.echoInfo;
    isarBeat               = rvIsarResult.isar.beat;
    isarBeatInfo           = rvIsarResult.isar.beatInfo;
    isarPost               = rvIsarResult.isar.post;
    isarSharedParameters   = rvIsarResult.isar.sharedParameters;
    if isfield(isarPost, 'motionDiagnostics')
        initialLosMotionDiagnostics = isarPost.motionDiagnostics;
    end
end
