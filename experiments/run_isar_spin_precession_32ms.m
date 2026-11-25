function results=run_isar_spin_precession_32ms(options)
%RUN_ISAR_SPIN_PRECESSION_32MS Exact echo + existing ISAR chain acceptance.
% Example: run_isar_spin_precession_32ms
% Options: Modes=["spin","coneSpin"], SpinFrequenciesHz=1:3,
% PrecessionFrequenciesHz=1:3, OutputDirectory=<timestamped results folder>.
% Resume=true reuses case MAT files in an explicitly selected output folder.
% TargetIndex=1:10 selects the corresponding original fullpath centre state;
% RecordTimeSec=1e-3 uses the original full record (default 200 us for target 1).
% PrecessionHalfAngleDeg=10, PulseCount=32, GeneratePlots=true are overridable.
% One representative 10.05 km cone uses the fullpath target-1 radial motion.
% Only its empty record tail is shortened: T=200 us, PRI=1 ms, N=32. All
% waveform, dechirp, RCMC, autofocus, and image settings remain the ISARpath
% defaults. No noise is added. The 5 cm criterion is worst point half-power
% full width, measured after replaying the FOUR-POINT scene's operators.
if nargin<1,options=struct();end
validateattributes(options,{'struct'},{'scalar'});
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'functions')); ensure_md_zt_paths();
modes=string(opt(options,'Modes',["spin","coneSpin"]));
for k=1:numel(modes)
    switch lower(modes(k))
        case "spin",modes(k)="spin";
        case "conespin",modes(k)="coneSpin";
        otherwise,error('RadarProject:ISARStudy:InvalidMode','Use spin or coneSpin.');
    end
end
spinFrequencies=opt(options,'SpinFrequenciesHz',1:3);
precessionFrequencies=opt(options,'PrecessionFrequenciesHz',1:3);
out=char(opt(options,'OutputDirectory',fullfile(root,'results', ...
    'isar_spin_precession_32ms',char(datetime('now','Format','yyyyMMdd_HHmmss')))));
if ~isfolder(out),mkdir(out);end
diary(fullfile(out,'run_log.txt')); diaryCleanup=onCleanup(@()diary('off')); %#ok<NASGU>
allRows=table();allPoints=table();cases=struct([]);
for mode=modes(:).'
    rates=precessionFrequencies;if mode=="spin",rates=0;end
    for fs=spinFrequencies(:).'
        for fp=rates(:).'
            fprintf('\nRunning %s: spin %.3g Hz, precession %.3g Hz, CPI %g ms\n', ...
                mode,fs,fp,opt(options,'PulseCount',32));
            label=sprintf('%s_fs%g_fp%g',mode,fs,fp);
            halfAngle=opt(options,'PrecessionHalfAngleDeg',10);
            pulseCount=opt(options,'PulseCount',32);
            if mode=="coneSpin"&&halfAngle~=10,label=sprintf('%s_a%g',label,halfAngle);end
            if pulseCount~=32,label=sprintf('%s_n%d',label,pulseCount);end
            caseFile=fullfile(out,[label '.mat']);
            if opt(options,'Resume',false)&&isfile(caseFile)
                load(caseFile,'entry','points','data');
                assert(data.cfg.T==opt(options,'RecordTimeSec',200e-6), ...
                    'Resume requires the same record length.');
                assert(data.cfg.Na==pulseCount&&data.cfg.precessionHalfAngleDeg==halfAngle, ...
                    'Resume requires the same aperture and half-cone angle.');
                if isfield(data,'inputMotion')
                    rangeList=[10050,20050,30050,40050,50050,65050,75050,85050,100050,120050];
                    assert(data.inputMotion.p0_I(1)==rangeList(opt(options,'TargetIndex',1)), ...
                        'Resume requires the same physical target.');
                else
                    assert(opt(options,'TargetIndex',1)==1,'Old cache represents target 1 only.');
                end
                fprintf('  Reusing completed case from this output directory.\n');
            else
                timer=tic;
                [entry,points,data]=run_case(root,mode,fs,fp,options);
                entry.runtimeSec=toc(timer);
                save(caseFile,'entry','points','data','-v7.3');
            end
            entry.targetIndex=opt(options,'TargetIndex',1);
            entry.recordTimeSec=data.cfg.T;
            allRows=[allRows;struct2table(entry)]; %#ok<AGROW>
            allPoints=[allPoints;points]; %#ok<AGROW>
            item=struct('label',label,'sceneImage',data.post.image{1}, ...
                'rangeAxisM',data.post.rangeAxisM{1},'crossAxisM',data.post.crossAxisM{1}, ...
                'profiles',{data.profiles},'widthsM',points.crossRange3dBWidthM, ...
                'aliasRisk',entry.microDopplerAliasRisk);
            cases=[cases;item]; %#ok<AGROW>
            writetable(allRows,fullfile(out,'summary.csv'));
            writetable(allPoints,fullfile(out,'point_widths.csv'));
            fprintf('  widths [cm]: %s; worst %.3f cm; alias risk %d; replay error %.3g\n', ...
                mat2str(points.crossRange3dBWidthM.'*100,4), ...
                entry.worstPointWidthM*100,entry.microDopplerAliasRisk,entry.replayRelativeError);
            clear data;
        end
    end
end
results=struct('summary',allRows,'points',allPoints,'outputDirectory',out, ...
    'matlabVersion',version,'criterion',"worst point power -3 dB full width <= 0.05 m", ...
    'noiseAdded',false,'cases',cases);
save(fullfile(out,'results.mat'),'results','-v7.3');
if opt(options,'GeneratePlots',true),plot_results(cases,out);end
disp(allRows(:,{'mode','spinHz','precessionHz','worstPointWidthM', ...
    'allPointWidthsWithin5cm','microDopplerAliasRisk','unambiguous5cm'}));
fprintf('Results saved to %s\n',out);
end

function [entry,points,data]=run_case(root,mode,fs,fp,options)
c=3e8;
targetIndexToSimulate=opt(options,'TargetIndex',1);
validateattributes(targetIndexToSimulate,{'numeric'},{'scalar','integer','>=',1,'<=',10});
rangeList=[10050,20050,30050,40050,50050,65050,75050,85050,100050,120050];
accelerationList=[-150,-120,-110,-108,-106,-103,-102,-101,-101,-101];
recordTimeSec=opt(options,'RecordTimeSec',200e-6);
cfg_ISARpath=struct('fs',1e9,'Tr',100e-6,'T',200e-6,'Ta',1e-3, ...
    'fL',11.8e9,'fH',17.8e9,'Na',32,'fs_scope',200e6, ...
    'echoMode',"fastTimeExact",'propagationSpeed',c,'startTimeSec',0, ...
    'bodyRotationMode',mode,'spinFrequencyHz',fs, ...
    'precessionFrequencyHz',max(1,fp),'precessionHalfAngleDeg',10, ...
    'omegaBody_B',[0 2*pi*fs 0],'fastTimeChunkSize',4096);
cfg_ISARpath.T=recordTimeSec;
cfg_ISARpath.Na=opt(options,'PulseCount',32);
validateattributes(cfg_ISARpath.Na,{'numeric'},{'scalar','integer','>=',2});
cfg_ISARpath.precessionHalfAngleDeg=opt(options,'PrecessionHalfAngleDeg',10);
cfg_ISARpath.frequencyBiasPatternHz=20e6+25e6*(mod(targetIndexToSimulate,2)==0);
cfg_ISARpath.B=cfg_ISARpath.fH-cfg_ISARpath.fL;
cfg_ISARpath.kr=cfg_ISARpath.B/cfg_ISARpath.Tr;
cfg_ISARpath.fc=(cfg_ISARpath.fH+cfg_ISARpath.fL)/2;
cfg_ISARpath.lambda=c/cfg_ISARpath.fc;
motion=struct('type','accelerated','referenceTimeSec',0, ...
    'p0_I',[rangeList(targetIndexToSimulate) 0 0],'v0_I',[3468 0 0], ...
    'a0_I',[accelerationList(targetIndexToSimulate) 0 0], ...
    'attitude0YprRad',[-3*pi/4 0 0]);
spec=struct('Ntip',1,'radiusM',.3,'phase0Rad',pi/6, ...
    'axialSeparationM',1.5,'alpha0',[1;.3;.3;.3]);
targets=scatterermodel.createLegacyDiskTarget(spec,motion);
radarState=struct('position_I',[0 0 0],'velocity_I',[0 0 0], ...
    'carrierFrequencyHz',cfg_ISARpath.fc,'wavelengthM',cfg_ISARpath.lambda, ...
    'dopplerSign',-1,'polarization','VV');
R_isar=rangeList(targetIndexToSimulate);v_isar=3468;
ISAR_PLOT_ENABLED=false;ISAR_QUALITY_ENABLED=true;
ISAR_QUALITY_EXPECTED_PEAK_COUNT=4;ISAR_QUALITY_PRINT_ENABLED=false;
run(which('multiuav_isarpath.m'));
% The total must be a coherent sum; measure each point using operators fit
% once to that sum. Re-estimating alignment per point is expressly avoided.
replayedSum=zeros(size(isarPost.phaseCompensatedHrrf{1}));
pointRecords=struct([]);profiles=cell(4,1);pointHrrf=cell(4,1);
for j=1:4
    request=isarBeatRequest;request.referenceStartTimeSec=referenceStartTimeSecISAR(1);
    pointBeat=radar.lfm_local_reference_beat( ...
        {isarEcho.perTarget(1).groups(j).signal},request,isarBeatOptions);
    response=replay_isar_point_response(pointBeat.total,cfg_ISARpath,isarPost,1,v_isar);
    replayedSum=replayedSum+response.hrrf;
    [width,profile]=point_width(response.hrrf,isarPost,1);
    profiles{j}=profile;
    pointHrrf{j}=response.hrrf;
    record=struct('mode',mode,'spinHz',fs,'precessionHz',fp,'pointIndex',j, ...
        'pointName',targets.groupNames(j),'crossRange3dBWidthM',width, ...
        'within5cm',isfinite(width)&&width<=.05,'mainlobeEnergyFraction',profile.mainlobeEnergyFraction);
    pointRecords=[pointRecords;record]; %#ok<AGROW>
end
points=struct2table(pointRecords);
replayError=norm(replayedSum-isarPost.phaseCompensatedHrrf{1},'fro')/ ...
    norm(isarPost.phaseCompensatedHrrf{1},'fro');
assert(replayError<1e-8,'Scene-to-point operator replay is not linear.');
% Dense exact geometry catches fractional-CPI Doppler extrema, including
% the intrapulse receive-time offset. Use fH for conservative sampling risk.
g=isarPost.crossRangeGeometry(1);
query=linspace(2*R_isar/cfg_ISARpath.propagationSpeed, ...
    (cfg_ISARpath.Na-1)*cfg_ISARpath.Ta+ ...
    2*(R_isar+v_isar*(cfg_ISARpath.Na-1)*cfg_ISARpath.Ta)/cfg_ISARpath.propagationSpeed+ ...
    cfg_ISARpath.Tr,257).';
state=scatterermodel.evaluateTarget(targetsISAR,radarStateISAR,query);
centre=state.centerPosition_I; centreV=state.centerVelocity_I;
centreRate=sum(centreV.*centre,2)./vecnorm(centre,2,2);
microDoppler=2*(state.rangeRateMps-centreRate.')/(c/cfg_ISARpath.fH);
span=max(microDoppler,[],1)-min(microDoppler,[],1);
aliasRisk=max(abs(microDoppler),[],'all')>=1/(2*cfg_ISARpath.Ta);
entry=struct('mode',mode,'spinHz',fs,'precessionHz',fp, ...
    'halfConeDeg',cfg_ISARpath.precessionHalfAngleDeg*(mode=="coneSpin"), ...
    'pulseCount',cfg_ISARpath.Na,'cpiSec',cfg_ISARpath.Na*cfg_ISARpath.Ta, ...
    'firstLastPulseSpanSec',(cfg_ISARpath.Na-1)*cfg_ISARpath.Ta,'prfHz',1/cfg_ISARpath.Ta, ...
    'effectiveRateRadps',g.effectiveRateRadps,'idealRayleighM',g.nominalResolutionM, ...
    'idealPower3dBWidthM',.885893*g.nominalResolutionM, ...
    'planeDriftDeg',g.maximumPlaneDriftDeg,'relativeRateVariation',g.relativeRateVariation, ...
    'meanPointWidthM',mean(points.crossRange3dBWidthM), ...
    'worstPointWidthM',max(points.crossRange3dBWidthM), ...
    'allPointWidthsWithin5cm',all(points.within5cm), ...
    'microDopplerMinHz',min(microDoppler,[],'all'), ...
    'microDopplerMaxHz',max(microDoppler,[],'all'), ...
    'maximumInstantaneousDopplerSpanHz',max(span), ...
    'microDopplerAliasRisk',aliasRisk,'unambiguous5cm',all(points.within5cm)&&~aliasRisk, ...
    'replayRelativeError',replayError);
data=struct('cfg',cfg_ISARpath,'post',isarPost,'quality',isarQuality, ...
    'inputMotion',motion, ...
    'profiles',{profiles},'pointHrrf',{pointHrrf}, ...
    'microDopplerHz',microDoppler,'microDopplerTimeSec',query);
end

function [width,profile]=point_width(H,post,k)
% 64x slow-time interpolation refines a width; it adds no aperture data.
A=abs(fftshift(fft(H,post.parameters.slowTimeZeroPaddingFactor*size(H,1),1),1));
[~,index]=max(A(:));[~,column]=ind2sub(size(A),index);
nfft=64*size(H,1);line=abs(fftshift(fft(H(:,column),nfft)));
[~,peak]=max(line);line=circshift(line,nfft/2+1-peak);line=line/max(line);
step=post.crossRangeGeometry(k).nominalResolutionM/64;
axisM=((0:nfft-1)-nfft/2)*step;
m=signal.resolution_measure_3db(line,axisM,'linear',struct(),false);
width=m.width;power=line.^2;
inside=axisM>=m.x_left & axisM<=m.x_right;
profile=struct('axisRelativeM',axisM,'amplitude',line, ...
    'mainlobeEnergyFraction',sum(power(inside))/sum(power));
end

function plot_results(cases,out)
n=numel(cases);cols=min(3,n);rows=ceil(n/cols);
fig=figure('Visible','off','Color','w','Position',[100 100 1300 300*rows]);
tiledlayout(rows,cols,'TileSpacing','compact','Padding','compact');
sgtitle('32 ms CPI, noiseless full chain; * indicates micro-Doppler alias risk');
for k=1:n
    nexttile;imagesc(cases(k).rangeAxisM,cases(k).crossAxisM, ...
        20*log10(max(cases(k).sceneImage,1e-2)));axis xy;clim([-40 0]);
    support=max(cases(k).sceneImage,[],1)>.05;
    rangeSupport=cases(k).rangeAxisM(support);
    if ~isempty(rangeSupport),xlim([min(rangeSupport)-.1,max(rangeSupport)+.1]);end
    xlabel('Range offset (m)');ylabel('Local cross-range (m)');
    title(case_title(cases(k)));
end
colormap(fig,'turbo');exportgraphics(fig,fullfile(out,'scene_images.png'),'Resolution',150);close(fig);
fig=figure('Visible','off','Color','w','Position',[100 100 1300 300*rows]);
tiledlayout(rows,cols,'TileSpacing','compact','Padding','compact');
sgtitle('Four-point scene operators applied to each point; * indicates alias risk');
for k=1:n
    nexttile;hold on;
    for j=1:4,p=cases(k).profiles{j};plot(p.axisRelativeM*100,20*log10(max(p.amplitude,1e-3)));end
    yline(-3,'k:');xline(-2.5,'k:');xline(2.5,'k:');grid on;
    xlim([-15 15]);ylim([-30 0]);xlabel('Offset from each peak (cm)');ylabel('Magnitude (dB)');
    title(sprintf('%s | worst %.2f cm',case_title(cases(k)),max(cases(k).widthsM)*100));
end
legend('Tip','Base 1','Base 2','Base 3','Location','southwest');
exportgraphics(fig,fullfile(out,'point_response_profiles.png'),'Resolution',150);close(fig);
valid=find(arrayfun(@(x) ~x.aliasRisk&&all(x.widthsM<=.05),cases));
if ~isempty(valid)
    [~,best]=min(arrayfun(@(x) max(x.widthsM),cases(valid)));k=valid(best);
    fig=figure('Visible','off','Color','w','Position',[100 100 1200 460]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile;imagesc(cases(k).rangeAxisM,cases(k).crossAxisM, ...
        20*log10(max(cases(k).sceneImage,1e-2)));axis xy;clim([-40 0]);
    support=max(cases(k).sceneImage,[],1)>.05;rangeSupport=cases(k).rangeAxisM(support);
    xlim([min(rangeSupport)-.1,max(rangeSupport)+.1]);
    xlabel('Range offset (m)');ylabel('Local cross-range (m)');
    title(case_title(cases(k)));colormap(fig,'turbo');colorbar;
    nexttile;hold on;
    names=["Tip","Base 1","Base 2","Base 3"];
    for j=1:4
        p=cases(k).profiles{j};plot(p.axisRelativeM*100,20*log10(max(p.amplitude,1e-3)), ...
            'LineWidth',1.3,'DisplayName',sprintf('%s: %.2f cm',names(j),cases(k).widthsM(j)*100));
    end
    yline(-3,'k:','HandleVisibility','off');grid on;xlim([-10 10]);ylim([-25 0]);
    xlabel('Offset from each peak (cm)');ylabel('Magnitude (dB)');
    title('32 ms CPI: measured point responses');legend('Location','southwest');
    exportgraphics(fig,fullfile(out,'recommended_case.png'),'Resolution',160);close(fig);
end
end

function titleText=case_title(item)
titleText=strrep(item.label,'_',' ');
if item.aliasRisk,titleText=[titleText ' *'];end
end

function value=opt(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
