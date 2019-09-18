load('Results/Lundstrom/Lundstrom_phaseInfo.mat');
load('Results/Lundstrom/Lundstrom_timeInfo.mat');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat', 'spkHistLengths');
load('Results/Lundstrom/alphas.mat');

spkHistLengths(1) = 10;

addpath ~/FigureComposer/
addpath ~/FigureComposer/matlab/


saveFigs = false;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 3;

figDir = 'Figs/FractionalDiff/';

%%
lw = 1.5;
lw2 = 1;

T = 4;
H = 26;

NR = 3;
NC = 3;

Hs = [8 length(spkHistLengths)];

colors = [0.6 0.6 0.6;
          0   0   0];


D = 4;

      
figure(2);
clf



% plot tau
subplot(NR,NC,1);
loglog(StimPeriods,nanmean(fits_time.tau_3AHP(:,D,1),2),'^-','LineWidth',lw,'Color',[0    0    1]);
hold on
plot(StimPeriods,nanmean(fits_time.tau_3AHP(:,D,2),2),'v-','LineWidth',lw,'Color',[0    0    1]);
for hh = 1:length(Hs)
    plot(StimPeriods,nanmean(fits_time.tau_glm_3AHP(:,D,Hs(hh),1),2),'^-','LineWidth',lw,'Color',colors(hh,:));
    plot(StimPeriods,nanmean(fits_time.tau_glm_3AHP(:,D,Hs(hh),2),2),'v-','LineWidth',lw,'Color',colors(hh,:));
end
xlabel('stimulus period (s)','FontSize',fontSizeLabel);
ylabel('time constant (s)','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
title('time constant','FontSize',fontSizeTitle);
hold off

% plot gain
subplot(NR,NC,2);
semilogx(StimPeriods,nanmean(fits_phase.gs_3AHP(:,D),2),'o-','LineWidth',lw,'Color',[0    0    1]);
hold on
XX = [log(StimPeriods(:)) ones(length(StimPeriods),1)];
b_hh = (XX'*XX)\(XX'*fits_phase.gs_3AHP(:,D));
Y_est = XX*b_hh;
plot(StimPeriods,nanmean(Y_est,2),'--','LineWidth',lw2,'Color',[0    0    1]);
xlabel('stimulus period (s)','FontSize',fontSizeLabel);
ylabel('gain (spks/s)','FontSize',fontSizeLabel);
title('gain','FontSize',fontSizeTitle);


for hh = 1:length(Hs)
    plot(StimPeriods,nanmean(fits_phase.gs_glm_3AHP(:,D,Hs(hh)),2),'o-','LineWidth',lw,'Color',colors(hh,:));
    b_glm = (XX'*XX)\(XX'*fits_phase.gs_glm_3AHP(:,D,Hs(hh)));
    Y_est = XX*b_glm;
    plot(StimPeriods,nanmean(Y_est,2),'--','LineWidth',lw2,'Color',colors(hh,:));
end
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off

% plot phase lead
subplot(NR,NC,3);
semilogx(StimPeriods,nanmean(rad2deg(fits_phase.phase_3AHP(:,D)-fits_phase.phase_Gain(:,D)),2),'o-','LineWidth',lw,'Color',[0    0    1]);
hold on
plot(StimPeriods([1 end]),nanmean(nanmean(rad2deg(fits_phase.phase_3AHP(:,D)-fits_phase.phase_Gain(:,D))),2)*[1 1],'--','LineWidth',lw2,'Color',[0    0    1]);
for hh = 1:length(Hs)
    plot(StimPeriods,nanmean(rad2deg(fits_phase.phase_glm_3AHP(:,D,Hs(hh))-fits_phase.phase_Gain(:,D)),2),'o-','LineWidth',lw,'Color',colors(hh,:));
    plot(StimPeriods([1 end]),nanmean(nanmean(rad2deg(fits_phase.phase_glm_3AHP(:,D,Hs(hh))-fits_phase.phase_Gain(:,D))),2)*[1 1],'--','LineWidth',lw2,'Color',colors(hh,:));
end
xlabel('stimulus period (s)','FontSize',fontSizeLabel);
ylabel('phase lead (deg)','FontSize',fontSizeLabel);
title('phase lead','FontSize',fontSizeTitle);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);

hold off

dd = 1e3;
% plot alphas: tau
subplot(NR,NC,4);

semilogx(spkHistLengths./dd,ones(size(spkHistLengths))*nanmean(alphas.filter.hh(D-1),2),'-','LineWidth',lw,'Color',[0 0 1]);
hold on
semilogx(spkHistLengths./dd,nanmean(alphas.filter.flat(:,D-1),2)  ,'LineWidth',lw,'Color',[0.8500    0.3250    0.0980]);
semilogx(spkHistLengths./dd,nanmean(alphas.filter.square(:,D-1),2)  ,'LineWidth',lw,'Color',[0    0 0]);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
xlim([10 16e3]./dd);
if(dd == 1e3)
    xlabel('spike history length (s)','FontSize',fontSizeLabel);
else
    xlabel('spike history length (ms)','FontSize',fontSizeLabel);
end
ylabel('estimated alpha','FontSize',fontSizeLabel);
ylim([ 0 0.22]);
hold off


% plot alphas: gain
subplot(NR,NC,5);


semilogx(spkHistLengths./dd,ones(size(spkHistLengths))*nanmean(alphas.gain.hh(D-1),2),'-','LineWidth',lw,'Color',[0 0 1]);
hold on
semilogx(spkHistLengths./dd,nanmean(alphas.gain.flat(:,D-1),2)  ,'LineWidth',lw,'Color',[0.8500    0.3250    0.0980]);
semilogx(spkHistLengths./dd,nanmean(alphas.gain.sine(:,D-1),2)  ,'LineWidth',lw,'Color',[0    0 0]);
xlim([10 16e3]./dd);
if(dd == 1e3)
    xlabel('spike history length (s)','FontSize',fontSizeLabel);
else
    xlabel('spike history length (ms)','FontSize',fontSizeLabel);
end
ylim([ 0 0.22]);
ylabel('estimated alpha','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off

% plot alphas: phase lead
subplot(NR,NC,6);

semilogx(spkHistLengths./dd,ones(size(spkHistLengths))*nanmean(alphas.phase.hh(D-1),2),'-','LineWidth',lw,'Color',[0 0 1]);
hold on
semilogx(spkHistLengths./dd,nanmean(alphas.phase.flat(:,D-1),2)  ,'LineWidth',lw,'Color',[0.8500    0.3250    0.0980]);
semilogx(spkHistLengths./dd,nanmean(alphas.phase.sine(:,D-1),2)  ,'LineWidth',lw,'Color',[0    0 0]);
xlim([10 16e3]./dd);
if(dd == 1e3)
    xlabel('spike history length (s)','FontSize',fontSizeLabel);
else
    xlabel('spike history length (ms)','FontSize',fontSizeLabel);
end
ylim([ 0 0.22]);
ylabel('estimated alpha','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off


% plot alphas by sigma: timescale
subplot(NR,NC,7);
hold on
plot(StimLevels(2:end),alphas.filter.hh,'o-','LineWidth',lw,'Color',[0 0 1]);
plot(StimLevels(2:end),alphas.filter.flat(end,:),'o-','LineWidth',lw,'Color',[0.8500    0.3250    0.0980]);
plot(StimLevels(2:end),alphas.filter.square(end,:),'o-','LineWidth',lw,'Color',[0 0 0]);
xlabel('sigma','FontSize',fontSizeLabel);
ylabel('time constant (s)','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
ylim([0.14 0.23]);
hold off

% plot alphas by sigma: gain
subplot(NR,NC,8);
hold on
plot(StimLevels(2:end),alphas.gain.hh,'o-','LineWidth',lw,'Color',[0 0 1]);
plot(StimLevels(2:end),alphas.gain.flat(end,:),'o-','LineWidth',lw,'Color',[0.8500    0.3250    0.0980]);
plot(StimLevels(2:end),alphas.gain.sine(end,:),'o-','LineWidth',lw,'Color',[0 0 0]);
xlabel('sigma','FontSize',fontSizeLabel);
ylabel('gain (spks/s)','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
ylim([0.14 0.23]);
hold off

% plot alphas by sigma: phase lead
subplot(NR,NC,9);
hold on
plot(StimLevels(2:end),alphas.phase.hh,'o-','LineWidth',lw,'Color',[0 0 1]);
plot(StimLevels(2:end),alphas.phase.flat(end,:),'o-','LineWidth',lw,'Color',[0.8500    0.3250    0.0980]);
plot(StimLevels(2:end),alphas.phase.sine(end,:),'o-','LineWidth',lw,'Color',[0 0 0]);
xlabel('sigma','FontSize',fontSizeLabel);
ylabel('phase lead (deg)','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
ylim([0.14 0.23]);
hold off



set(gcf,'PaperUnits','inches','PaperSize',[NC NR]*PlotSize,'PaperPosition',[0 0 NC NR]*PlotSize);
if(saveFigs)
    matfig2fyp(gcf,sprintf('%s/fig2_alphas_raw.fyp', figDir));
end


