function [name,audit,fit]=isar_tip_colour_overlay(one,out,rho)
%ISAR_TIP_COLOUR_OVERLAY Merge the measured spectrum, theory and peak fit.
% Fit ONLY the saved extracted peaks. The theoretical trajectory is never
% an input to the fit, and the original peak/theory correlation is retained.
q=one.scene;X=one.tracks;G=one.truth;
t=X.receiveMidpointSec;observed=X.extractedTipPeakHz;
take=t<=3;dense=G.slowTimeSec<=3;grid=(0:.001:3)';
step=median(diff(t));nyquist=1/(2*step);frequencies=.5:.5:4.5;
assert(numel(t)==50&&nnz(take)==30&&max(frequencies)<nyquist);
origin=mean(t);A=fit_matrix(t,origin,frequencies);
coefficient=A\observed;
fitAtPeaks=A*coefficient;curve=fit_matrix(grid,origin,frequencies)*coefficient;
fit=struct('caseId',q.caseId,'method',"HARMONIC_LS_FROM_EXTRACTED_PEAKS", ...
    'frequenciesHz',frequencies,'timeOriginSec',origin,'coefficient',coefficient, ...
    'fitCpiCount',numel(t),'fitUsesTheoreticalTrajectory',false, ...
    'affectsProcessingOrOriginalCorrelation',false);
residual=fitAtPeaks-observed;
name=string(fullfile(out,q.caseId+"_colour_3s.png"));
audit=table(q.caseId,rho,numel(t),nnz(take),size(A,2),nyquist, ...
    sqrt(mean(residual(take).^2)),sqrt(mean(residual.^2)), ...
    max(abs(residual(take))),name, ...
    'VariableNames',{'caseId','originalPeakTheoryR3s','fitCpiCount', ...
    'displayedCpiCount','fitParameterCount','peakSamplingNyquistHz', ...
    'fitRms3sHz','fitRms5sHz','maximumFitResidual3sHz','colourImage'});

colours=turbo(256);
span=max(50,25*ceil((max(abs([G.tipPhysicalDopplerHz;observed]))+25)/25));
span=min(span,250);tick=25;if span>100,tick=50;end
top=sprintf('自旋 %g Hz，进动 %g Hz，半锥角 %g°', ...
    q.spinFrequencyHz,q.precessionFrequencyHz,q.precessionHalfAngleDeg);
if q.swingAmplitudeDeg>0
    top={top,sprintf('固定Z轴摆动 ±%g°，%g Hz',q.swingAmplitudeDeg,q.swingFrequencyHz)};
end
fig=figure('Visible','off','Color','w','Position',[30 30 1440 1100]);
guard=onCleanup(@()close(fig)); %#ok<NASGU>
layout=tiledlayout(fig,1,1,'Padding','compact');ax=nexttile(layout);hold(ax,'on');
power=one.tipGatePower(:,take);
db=10*log10(max(power./max(power,[],1),1e-12));
im=imagesc(ax,t(take),one.dopplerHz,db);im.AlphaData=1;im.Tag='saved_tip_spectrum';
ax.Color=colours(1,:);axis(ax,'xy');colormap(ax,colours);clim(ax,[-10 0]);
h1=plot(ax,G.slowTimeSec(dense),G.tipPhysicalDopplerHz(dense),'w--', ...
    'LineWidth',2.2,'Tag','physical_theory');
h2=plot(ax,grid,curve,'-','Color',[0 .88 1], ...
    'LineWidth',2.2,'Tag','peak_harmonic_fit');
h3=plot(ax,t(take),observed(take),'o','LineStyle','none', ...
    'MarkerSize',7.5,'MarkerEdgeColor',[.05 .05 .05], ...
    'MarkerFaceColor','w','LineWidth',1.3,'Tag','extracted_peak_markers');
bar=colorbar(ax);bar.FontSize=30;bar.Ticks=[-10 -5 0];
% Use a dark dashed legend sample so the white theory line remains legible
% in a transparent legend against the figure's white margin.
theoryKey=plot(ax,NaN,NaN,'--','Color',[.25 .25 .25], ...
    'LineWidth',h1.LineWidth,'Tag','theory_legend_sample');
legend(ax,[theoryKey h2 h3],{'理论轨迹','峰值拟合','提取峰值'}, ...
    'Location','southoutside','Orientation','horizontal','FontSize',30, ...
    'FontName','Microsoft YaHei','Color','none','Box','off', ...
    'TextColor','k','AutoUpdate','off');
xlim(ax,[0 3]);xticks(ax,0:.5:3);ylim(ax,[-span span]);yticks(ax,-span:tick:span);
title(ax,top,'FontSize',30,'FontWeight','normal','FontName','Microsoft YaHei');
subtitle(ax,sprintf('R = %.6f',rho),'FontSize',30,'FontName','Microsoft YaHei');
ax.FontSize=30;ax.FontName='Microsoft YaHei';ax.LineWidth=1.2;ax.TickDir='out';
ax.XLabel.String='';ax.YLabel.String='';set(findall(fig,'Type','text'),'FontSize',30);
fonts=arrayfun(@(h) double(h.FontSize),findall(fig,'-property','FontSize'));
assert(all(abs(fonts-30)<1e-9)&&isempty(ax.XLabel.String)&&isempty(ax.YLabel.String));
assert(numel(h3.XData)==30&&numel(h2.XData)==3001&&im.AlphaData==1);
drawnow;
savefig(fig,replace(name,'.png','.fig'));
exportgraphics(fig,name,'Resolution',120);
end

function A=fit_matrix(t,origin,frequencies)
t=t(:);phase=2*pi*t*frequencies;
A=[ones(size(t)),t-origin,cos(phase),sin(phase)];
end
