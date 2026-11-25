function test_isar_ca_run_contract(parent)
% Small synthetic data only. No RF, image formation or identification run.
if nargin == 0, parent = tempname; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'functions'));
t = (0:0.1:2).';
r = 1+2*t+1.5*t.^2;
v = 2+3*t;
cfg = struct('propagationSpeed',3e8,'kr',6e13);
g = isar_golden_target_config();
infos = cell(3,1);
modes = ["measured","exact","measured"];
for j=1:numel(modes)
    info = isar_ca_run_info(parent,'test',"tip",modes(j),0.02,0.3);
    infos{j} = info;
    rec = isar_cone_measurement_record(t,r,v,struct('Mode',modes(j),'Seed',17));
    a = isar_ca_record_acceleration(rec);
    if modes(j) == "exact"
        assert(abs(a-3)<1e-10 && info.appliedRangeSigmaM==0 && info.appliedVelocitySigmaMps==0);
    else
        assert(abs(a-3)>1e-8 && info.appliedRangeSigmaM==0.02 && info.appliedVelocitySigmaMps==0.3);
    end
    rec2 = rec;
    rec2.rangeM = rec2.rangeM+0.5*0.1*t.^2;
    assert(abs(isar_ca_record_acceleration(rec2)-a-0.1)<1e-10);
    assert(info.rangeModel=="tip" && info.velocityModel=="centre");
    % Audit every recorded offset on two distinct frames.
    for it=1:numel(g.beatOffsetHz)
        for k=[1 numel(t)]
            bias = g.beatOffsetHz(it);
            pr = isar_cone_processing_prior(rec,k,a,bias,cfg);
            row = isar_ca_lo_audit(g.targetId(it),k,bias,pr);
            assert(row.frameIndex==k && row.rangeM==rec.rangeM(k) && row.epochSec==t(k));
            assert(row.targetId==it-1 && row.usedBiasHz==bias && row.goldenBiasHz==bias);
            assert(abs(row.loStartSec-(2*row.rangeM/cfg.propagationSpeed-bias/cfg.kr))<1e-15);
        end
    end
    % A reporting error must remain visible, not become "expected == used".
    row = isar_ca_lo_audit(0,k,pr.beatOffsetHz+1,pr);
    assert(row.usedBiasHz==pr.beatOffsetHz && row.usedBiasHz~=row.goldenBiasHz);
    tbl = isar_ca_tag_summary(table(rec.seed,'VariableNames',{'recordSeed'}),info);
    writetable(tbl,fullfile(info.outputDirectory,'summary.csv'));
    save(fullfile(info.outputDirectory,'summary.mat'),'info','tbl');
end
for j=1:numel(infos)
    info = infos{j};
    data = load(fullfile(info.outputDirectory,'summary.mat'));
    csv = readtable(fullfile(info.outputDirectory,'summary.csv'),'TextType','string');
    assert(isequal(data.info,info) && csv.runId==info.runId && csv.recordMode==modes(j));
    assert(csv.recordSeed==17 && csv.velocityModel=="centre" && csv.rangeModel=="tip");
end
assert(numel(unique(cellfun(@(s) char(s.runId),infos,'UniformOutput',false)))==3);
fprintf('test_isar_ca_run_contract passed: %s\n',parent);
end
