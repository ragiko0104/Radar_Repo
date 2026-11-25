function result=run_isar_tip_motion_retest(options)
%RUN_ISAR_TIP_MOTION_RETEST Per-scene native-RF tip analysis and stage gate.
% Stage 2/3 cannot start unless stage 1 has finished and its tip-position
% sampling audit passes. A spin Euler factor is not a tip motion frequency.
if nargin<1,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));
[~,~,~,template]=isar_grid_setup();
out=char(opt(options,'OutputDirectory',fullfile(root,'results','isar_tip_motion_retest_20261020')));
if ~isfolder(out),mkdir(out);end
groups=opt(options,'Groups',[1 2]);assert(issorted(groups));
summary=table();manifest=table();status="COMPLETE_REQUESTED_GROUPS";
for group=groups
    gateFile=fullfile(out,'first_group_sampling_gate.csv');
    if group>1
        assert(isfile(gateFile),'Stage 1 must finish before later stages.');
        gate=readtable(gateFile,'TextType','string','Delimiter',',','VariableNamingRule','preserve');
        if ~all(gate.tipPositionNyquistSatisfied)
            status="STOPPED_AFTER_GROUP1";break;
        end
    end
    cases=make_cases(template(1,:),group,options);
    groupOut=fullfile(out,sprintf('group%d',group));
    if ~isfolder(groupOut),mkdir(groupOut);end
    writetable(cases,fullfile(groupOut,'physical_scenes.csv'));
    run_isar_tip_slowtime_tracks(struct('Cases',cases,'OutputDirectory',groupOut, ...
        'NativeWorkers',opt(options,'NativeWorkers',8),'DurationSec',5, ...
        'SavedDopplerHalfSpanHz',250,'GenerateFigures',false,'Resume',true, ...
        'CaseCompletedFcn',@(one) render_completed_scene(one,groupOut)));
    [stats,images]=finish_group(groupOut);
    summary=[summary;stats];manifest=[manifest;images]; %#ok<AGROW>
    if group==1
        gate=audit_first_group(cases,groupOut);
        writetable(gate,gateFile);
        if ~all(gate.tipPositionNyquistSatisfied)
            status="STOPPED_AFTER_GROUP1";
            fprintf('FIRST_GROUP_GATE FAILED: stop all later groups.\n');break;
        end
        fprintf('FIRST_GROUP_GATE PASSED: tip 3 Hz < 5 Hz display Nyquist; later groups may proceed.\n');
    end
end
% Include previously completed groups when a later invocation handles group 3.
allSummary=table();allManifest=table();
for group=1:3
    folder=fullfile(out,sprintf('group%d',group));
    if isfile(fullfile(folder,'scene_similarity.csv'))
        allSummary=[allSummary;readtable(fullfile(folder,'scene_similarity.csv'), ...
            'TextType','string','Delimiter',',','VariableNamingRule','preserve')]; %#ok<AGROW>
        allManifest=[allManifest;readtable(fullfile(folder,'scene_figures.csv'), ...
            'TextType','string','Delimiter',',','VariableNamingRule','preserve')]; %#ok<AGROW>
    end
end
writetable(allSummary,fullfile(out,'all_scene_similarity.csv'));
writetable(allManifest,fullfile(out,'all_scene_figures.csv'));
result=struct('status',status,'outputDirectory',string(out),'summary',allSummary, ...
    'figures',allManifest,'requestedGroups',groups,'rawEchoSaved',false);
save(fullfile(out,'batch_summary.mat'),'result');disp(allSummary);
fprintf('BATCH %s | %s\n',status,out);
end

function cases=make_cases(template,group,options)
cases=template([],:);
if group==1
    spin=3;prec=3;angles=[1 3];swing=0;
elseif group==2
    spin=[2 3.5];prec=[2 3.5];angles=[5 7];swing=0;
else
    assert(group==3&&isfield(options,'ThirdConingHalfAnglesDeg')&& ...
        isfield(options,'ThirdSwingAmplitudesDeg'),'Group 3 needs the clarified angle definitions.');
    spin=2;prec=3.5;angles=options.ThirdConingHalfAnglesDeg;swing=options.ThirdSwingAmplitudesDeg;
end
for fs=spin
for fp=prec
for alpha=angles
for amplitude=swing
    q=template;q.family=group;q.spinFrequencyHz=fs;q.precessionFrequencyHz=fp;
    q.precessionHalfAngleDeg=alpha;q.swingAmplitudeDeg=amplitude;
    q.swingFrequencyHz=double(group==3);q.mode="coneSpin";
    if group==3,q.mode="coneSwingSpin";end
    q.caseId=replace(string(sprintf('R%d_s%g_p%g_a%g_w%g',group,fs,fp,alpha,amplitude)),'.','p');
    cases=[cases;q]; %#ok<AGROW>
end
end
end
end
end

function [summary,manifest]=finish_group(out)
d=load(fullfile(out,'summary.mat'),'result');cases=d.result.cases;
summary=table();manifest=table();audit=table();
for ci=1:height(cases)
    q=cases(ci,:);data=load(fullfile(out,q.caseId+"_tip_result.mat"),'one');one=data.one;
    X=one.tracks;G=one.truth;take=X.receiveMidpointSec<=3;
    delta=interp1(G.slowTimeSec,G.tipRelativeCentreRangeM,X.receiveMidpointSec,'pchip');
    coupling=-4*6e13*(3400+100*X.receiveMidpointSec).*delta/(3e8)^2;
    theory=X.physicalTipMidpointHz;observed=X.extractedTipPeakHz;
    rho=corrcoef(theory(take),observed(take));rho5=corrcoef(theory,observed);
    corrected=observed-coupling;rhoC=corrcoef(theory(take),corrected(take));
    row=q;row.cpiCount=height(X);row.displayedCpiCount=nnz(take);
    row.zeroLagPearsonR3s=rho(1,2);row.zeroLagPearsonR5s=rho5(1,2);
    row.pearsonAfterCoupling3s=rhoC(1,2);
    row.meanBias3sHz=mean(observed(take)-theory(take));
    row.rawRmse3sHz=sqrt(mean((observed(take)-theory(take)).^2));
    row.rmseAfterCoupling3sHz=sqrt(mean((corrected(take)-theory(take)).^2));
    row.maximumTipPhysicalDopplerHz=max(abs(G.tipPhysicalDopplerHz));
    row.maximumTailLeakageDb=max(X.tailLeakageRelativeTipDb);
    row.maximumTailInducedPeakDifferenceHz=max(abs(X.tailInducedPeakDifferenceHz));
    row.tipMainPeakCoveredEveryCpi=all(X.tipMainPeakInGate);
    row.minimumTipTailRangeSeparationM=one.minimumTipTailRangeSeparationM;
    summary=[summary;row]; %#ok<AGROW>
    X.couplingPredictionHz=coupling;X.theoryPlusCouplingHz=theory+coupling;
    audit=[audit;X]; %#ok<AGROW>
    names=render_scene(one,out,rho(1,2));
    manifest=[manifest;table(q.caseId,q.family,names(1),names(2), ...
        'VariableNames',{'caseId','group','colourImage','comparisonImage'})]; %#ok<AGROW>
end
writetable(summary,fullfile(out,'scene_similarity.csv'));
writetable(manifest,fullfile(out,'scene_figures.csv'));
writetable(audit,fullfile(out,'scene_coupling_audit.csv'));
end

function names=render_scene(one,out,rho)
q=one.scene;X=one.tracks;G=one.truth;
names=[string(fullfile(out,q.caseId+"_colour_3s.png")), ...
    string(fullfile(out,q.caseId+"_comparison_3s.png"))];
if all(isfile(names)),return;end
take=X.receiveMidpointSec<=3;dense=G.slowTimeSec<=3;
colours=turbo(256);span=max(50,25*ceil((max(abs([G.tipPhysicalDopplerHz;X.extractedTipPeakHz]))+25)/25));
span=min(span,250);tick=25;if span>100,tick=50;end
top=sprintf('自旋 %g Hz，进动 %g Hz，半锥角 %g°', ...
    q.spinFrequencyHz,q.precessionFrequencyHz,q.precessionHalfAngleDeg);
if q.swingAmplitudeDeg>0
    top={top,sprintf('固定Z轴摆动 ±%g°，%g Hz',q.swingAmplitudeDeg,q.swingFrequencyHz)};
end
for kind=1:2
    if kind==1
        names(1)=isar_tip_colour_overlay(one,out,rho);
        continue;
    end
    fig=figure('Visible','off','Color','w','Position',[30 30 1440 1100]);
    layout=tiledlayout(fig,1,1,'Padding','compact');ax=nexttile(layout);hold(ax,'on');
        coupling=-4*6e13*(3400+100*G.slowTimeSec).*G.tipRelativeCentreRangeM/(3e8)^2;
        h1=plot(ax,G.slowTimeSec(dense),G.tipPhysicalDopplerHz(dense),'k-','LineWidth',2);
        h2=plot(ax,G.slowTimeSec(dense),G.tipPhysicalDopplerHz(dense)+coupling(dense), ...
            '--','Color',[.05 .35 .75],'LineWidth',2);
        h3=plot(ax,X.receiveMidpointSec(take),X.extractedTipPeakHz(take), ...
            'o','Color',[.82 .12 .12],'LineWidth',1.5,'MarkerSize',8);
        legend(ax,[h1 h2 h3],{'物理理论','理论＋耦合近似','处理峰值'}, ...
            'Location','southoutside','Orientation','horizontal','FontSize',30, ...
            'FontName','Microsoft YaHei');grid(ax,'on');
        suffix="comparison_3s";subtitleText=sprintf('零时延 Pearson r = %.6f（s，Hz）',rho);
    xlim(ax,[0 3]);xticks(ax,0:.5:3);ylim(ax,[-span span]);yticks(ax,-span:tick:span);
    title(ax,top,'FontSize',30,'FontWeight','normal','FontName','Microsoft YaHei');
    subtitle(ax,subtitleText,'FontSize',30,'FontName','Microsoft YaHei');
    ax.FontSize=30;ax.FontName='Microsoft YaHei';ax.LineWidth=1.2;ax.TickDir='out';
    ax.XLabel.String='';ax.YLabel.String='';set(findall(fig,'Type','text'),'FontSize',30);
    fonts=arrayfun(@(h) double(h.FontSize),findall(fig,'-property','FontSize'));
    assert(all(abs(fonts-30)<1e-9)&&isempty(ax.XLabel.String)&&isempty(ax.YLabel.String));
    stem=q.caseId+"_"+suffix;names(kind)=string(fullfile(out,stem+".png"));
    savefig(fig,fullfile(out,stem+".fig"));exportgraphics(fig,names(kind),'Resolution',120);close(fig);
end
end

function gate=audit_first_group(cases,out)
[cfg,sensor,base,~]=isar_grid_setup();gate=table();
t=(0:1999)'*.001;samplingHz=10;
for k=1:height(cases)
    q=cases(k,:);work=cfg;work.bodyRotationMode=q.mode;
    work.spinFrequencyHz=q.spinFrequencyHz;work.precessionFrequencyHz=q.precessionFrequencyHz;
    work.precessionHalfAngleDeg=q.precessionHalfAngleDeg;
    target=isar_apply_micro_motion(base,work);state=scatterermodel.evaluateTarget(target,sensor,t);
    tip=reshape(state.position_I(1,:,:),3,[]).'-state.centerPosition_I;
    tip=tip-mean(tip,1);amplitude=sqrt(sum(abs(fft(tip,[],1)).^2,2));
    amplitude(1)=0;[~,index]=max(amplitude(1:1001));dominant=(index-1)/2;
    % For no-swing P1=(0,-L/2,0), Ry(psi)*P1=P1 exactly; its Cartesian
    % motion consists solely of DC and the coning frequency, not fs+fp.
    maxTipFrequency=q.precessionFrequencyHz;
    assert(abs(dominant-maxTipFrequency)<1e-9);
    data=load(fullfile(out,q.caseId+"_tip_result.mat"),'one');
    X=data.one.tracks;
    gate=[gate;table(q.caseId,maxTipFrequency,dominant,samplingHz,samplingHz/2, ...
        samplingHz/maxTipFrequency,2*maxTipFrequency<samplingHz, ...
        max(abs(data.one.truth.tipPhysicalDopplerHz)), ...
        all(abs(X.extractedTipPeakHz)<500), ...
        'VariableNames',{'caseId','maximumTipPositionFrequencyHz','denseGeometryDominantHz', ...
        'imageSamplingHz','imageNyquistHz','samplesPerMotionCycle', ...
        'tipPositionNyquistSatisfied','maximumInstantaneousTipDopplerHz','tipDopplerWithinPrfNyquist'})]; %#ok<AGROW>
end
disp(gate);
end

function render_completed_scene(one,out)
take=one.tracks.receiveMidpointSec<=3;
rho=corrcoef(one.tracks.physicalTipMidpointHz(take),one.tracks.extractedTipPeakHz(take));
render_scene(one,out,rho(1,2));
fprintf('SCENE FIGURES READY %s | r3s %.9f\n',one.scene.caseId,rho(1,2));
end

function value=opt(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
