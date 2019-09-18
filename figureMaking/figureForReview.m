load Results/Lundstrom/Lundstrom_phaseInfo.mat
load('/home/latimerk/gitCode/gainscaling/Results/Lundstrom/alphas.mat', 'spkHistLengths','StimPeriods')
spkHistLengths_lundstrom = spkHistLengths;
load('/home/latimerk/gitCode/gainscaling/Results/Mease/Mease_gainInfo.mat')

fontSizeTitle = 12;
fontSizeLabel = 10;
fontSizeAxis  = 10;
fontSizeLegend = 10;

figure(1);
clf



p_b = histc(randn(2e6,1),fits_gain.bins);
vv = p_b >= 1e3;
p_b = p_b./sum(p_b);

%% fr curves
NN = 15;
KK = 15;
TT = 5;
HH = 10;

colors = linspace(0.8,0,4)'*[1 1 1];

subplot(1,3,1);
for ii = 1:4
    xx = fits_gain.bins*StimLevels(ii);
    
    fr = fits_gain.glm.p_spk(:,NN,KK,ii,TT,HH)./p_b;
    fr(~vv) = nan;
    
    semilogy(xx,fr,'Color',colors(ii,:),'LineWidth',1.5);
    hold on
end

xlabel('$$s$$','FontSize',fontSizeLabel,'Interpreter','Latex')
ylabel('input-output $$\hat{R}_\sigma[s]$$','FontSize',fontSizeLabel,'Interpreter','Latex')

set(gca,'TickDir','out','fontSize',fontSizeAxis,'box','off');
hold off


%% scaled fr curves
subplot(1,3,2);
for ii = 1:4
    xx = fits_gain.bins;
    
    fr = fits_gain.glm.p_spk(:,NN,KK,ii,TT,HH)./p_b;
    fr(~vv) = nan;
    
    semilogy(xx,fr,'Color',colors(ii,:),'LineWidth',1.5);
    hold on
end
xlabel('$$\hat{s}$$','FontSize',fontSizeLabel,'Interpreter','Latex')
ylabel('scaled input-output $$\hat{R}_\sigma[\hat{s}]$$','FontSize',fontSizeLabel,'Interpreter','Latex')

set(gca,'TickDir','out','fontSize',fontSizeAxis,'box','off');
hold off


%% phase lead
subplot(2,3,3);
LL = 4;
TT2 = 4;

semilogx(StimPeriods([1 end]),mean(rad2deg(fits_phase.phase_3AHP(:,LL)-fits_phase.phase_Gain(:,LL)))*[1 1],'--','Color',[0 0 0],'LineWidth',1);
hold on
semilogx(StimPeriods,rad2deg(fits_phase.phase_3AHP(:,LL)-fits_phase.phase_Gain(:,LL)),'o-','Color',[0 0 0],'LineWidth',1.5);

semilogx(StimPeriods([1 end]),mean(rad2deg(fits_phase.phase_glm_3AHP(:,LL,HH2)-fits_phase.phase_Gain(:,LL)))*[1 1],'--','Color',[0 0.4470 0.7410],'LineWidth',1);
semilogx(StimPeriods,rad2deg(fits_phase.phase_glm_3AHP(:,LL,HH2)-fits_phase.phase_Gain(:,LL)),'o-','Color',[0 0.4470 0.7410],'LineWidth',1.5);

ylabel('phase lead (\circ)','fontSize',fontSizeLabel)

hold on

set(gca,'TickDir','out','fontSize',fontSizeAxis,'box','off');
hold off

%% phase gain
subplot(2,3,6);

XX = [log(StimPeriods(:)) ones(length(StimPeriods),1)];
bb_hh = (XX'*XX)\(XX'*log(fits_phase.gs_3AHP(:,LL)));
bb_glm = (XX'*XX)\(XX'*log(fits_phase.gs_glm_3AHP(:,LL,HH2)));
Y_hat_hh = exp(XX*bb_hh);
Y_hat_glm = exp(XX*bb_glm);


loglog(StimPeriods,Y_hat_hh,'--','Color',[0 0 0],'LineWidth',1);
hold on
loglog(StimPeriods,fits_phase.gs_3AHP(:,LL),'o-','Color',[0 0 0],'LineWidth',1.5);

loglog(StimPeriods,Y_hat_glm,'--','Color',[0 0.4470 0.7410],'LineWidth',1);
loglog(StimPeriods,fits_phase.gs_glm_3AHP(:,LL,HH2),'o-','Color',[0 0.4470 0.7410],'LineWidth',1.5);

ylabel('gain (spks/S)','fontSize',fontSizeLabel)
xlabel('stimulus period (S)','fontSize',fontSizeLabel)

hold on

set(gca,'TickDir','out','fontSize',fontSizeAxis,'box','off');
hold off

%%

set(gcf,'PaperSize',[9 3],'PaperPosition',[0 0 9 3]);

addpath ~/gitPublic/FigureComposer/matlab/
saveas(gcf,'Figs/figureForReview.fig');
savefigasis(gcf,'~/gitCode/gainscaling/Figs/figureForReview2.fig');