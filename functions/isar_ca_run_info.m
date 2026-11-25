function info = isar_ca_run_info(parent, entry, rangeModel, mode, sigmaR, sigmaV)
% OutputDirectory is a parent directory. Never reuse a previous run.
entry = string(entry);
rangeModel = string(rangeModel);
mode = string(mode);
if ~isscalar(mode) || ~ismember(mode, ["exact","measured"]) || ...
        ~isscalar(rangeModel) || ~ismember(rangeModel, ["centre","tip"])
    error('ISAR:Chain:RunMode', 'Invalid measurement mode or range source.');
end
if ~isfolder(parent)
    [ok,msg] = mkdir(parent);
    if ~ok, error('ISAR:Chain:RunDirectory','%s',msg); end
end
[~,nonce] = fileparts(tempname(parent));
runId = strjoin([entry,rangeModel,mode,string(nonce)],"_");
outDir = fullfile(parent,char(runId));
[ok,msg,id] = mkdir(outDir);
if ~ok || ~isempty(id)
    error('ISAR:Chain:RunDirectory','Could not reserve fresh directory %s: %s',outDir,msg);
end
% Keep model identity separate from whether noise was actually applied.
info = struct('runId',runId,'outputDirectory',string(outDir),'entry',entry, ...
    'rangeModel',rangeModel,'velocityModel',"centre",'recordMode',mode, ...
    'accelerationSource',"fit:record.rangeM",'configuredRangeSigmaM',sigmaR, ...
    'configuredVelocitySigmaMps',sigmaV,'appliedRangeSigmaM',sigmaR*(mode=="measured"), ...
    'appliedVelocitySigmaMps',sigmaV*(mode=="measured"));
end
