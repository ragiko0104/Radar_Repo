%MULTIUAV_ISARPATH ISAR 成像链路纯计算兼容薄脚本（调用 isarpath_run）。
% 推荐新代码直接调用 rv_isar_run(simCfg) 或 isarpath_run(scene, priorMotion)。
% 绘图请调用 code/plot_rv_isar_results.m，评估请调用 code/assess_rv_isar_results.m。

isarTruth__ = struct('R0', R_isar(:), 'v_radial0', v_isar(:));
if exist('initialLosMotionDiagnostics', 'var') && isstruct(initialLosMotionDiagnostics)
    isarTruth__.motionDiagnostics = initialLosMotionDiagnostics;
end

isarSceneInput__ = struct( ...
    'cfgISAR',    cfg_ISARpath, ...
    'targets',    targets, ...
    'radarState', radarState, ...
    'truth',      isarTruth__);

isarPrior__ = struct('rangeM', R_isar(:), 'velocityMps', v_isar(:));
isarRes__   = isarpath_run(isarSceneInput__, isarPrior__);

targetsISAR            = isarRes__.targetsISAR;
isarMotion             = isarRes__.motion;
isarCrossRangeGeometry = isarRes__.crossRangeGeometry;
isarScattererState     = isarRes__.scattererState;
isarEcho               = isarRes__.echo;
isarEchoInfo           = isarRes__.echoInfo;
isarBeat               = isarRes__.beat;
isarBeatInfo           = isarRes__.beatInfo;
isarPost               = isarRes__.post;
isarSharedParameters   = isarRes__.sharedParameters;
R_prior_ISAR           = isarRes__.scene.priorRangeM;
v_prior_ISAR           = isarRes__.scene.priorVelocityMps;
pulseTimeSecISAR       = isarRes__.scene.pulseTimeSec;
fastTimeSecISAR        = isarRes__.scene.fastTimeSec;
radarStateISAR         = isarRes__.scene.radarState;
if isfield(isarPost, 'motionDiagnostics')
    initialLosMotionDiagnostics = isarPost.motionDiagnostics;
end

clear isarTruth__ isarSceneInput__ isarPrior__ isarRes__
