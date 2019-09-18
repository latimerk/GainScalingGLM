load Results/Mease/Mease_gainInfo.mat
load Results/Mease/MeaseGLMs.mat
load('/home/latimerk/gitCode/GainScaling/Results/Mease/gainDistance.mat')
load('/home/latimerk/gitCode/GainScaling/Results/Mease/glmFitMetrics.mat')
figDir = 'Figs/GainScaling/';

addpath ~/FigureComposer/
addpath ~/FigureComposer/matlab/

addpath figureMaking/

saveFigs = false;


fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 3;
         
spkHistLengths(1) = 10;

%%
T = 5;
T2 = 1;
D = 3;
Hs = [1 3 5 16];
NH = length(Hs);
lims = [600 2000]+[-1 1]*60;

cols = repmat(linspace(0.8,0,NH)',[1 3]);
sta_cs = [0    0    0;
          80   0  5;
          160  0  10;
          240 0 15]./255;
            


Gs = 600:100:2000;

[G_Ks,G_Nas] = meshgrid(Gs,Gs);

ratios  = (Gs./Gs')';
ratios(isnan(gainDistance.hh.EM(:,:,1))) = nan;
[ratios,order] = sort(ratios(:));

G_Nas = G_Nas(order);
G_Ks = G_Ks(order);

NR = 2;
NC = 3;


gg_hh = gainDistance.hh.EM(:,:,D);
gg_hh = gg_hh(order);

legendObjs = cell(NH+1,1);
legendObjs{1} = 'Hodgkin-Huxley';

gg_glm = zeros(length(gg_hh),NH);
gg_glm2 = zeros(length(gg_hh),NH);
for ii = 1:NH
    gg_c = gainDistance.glm.EM(:,:,D,Hs(ii),T);
    gg_glm(:,ii) = gg_c(order);
    gg_c = gainDistance.glm.EM(:,:,D,Hs(ii),T2);
    gg_glm2(:,ii) = gg_c(order);
    legendObjs{ii+1} = sprintf('h_{spk} = %d ms',spkHistLengths(Hs(ii)));
end

figure(1);
clf;

subplot(NR,NC,2);
hold on
plot(ratios,gg_hh,'linewidth',2,'color',[0 0 1]);

for ii = 1:NH
    plot(ratios,gg_glm(:,ii),'linewidth',1,'color',cols(ii,:));
end

legend(legendObjs,'FontSize',fontSizeLabel);
xlabel('G_na/G_k');
ylabel('gain scaling:','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off


subplot(NR,NC,2+NC);
hold on
plot(ratios,gg_hh,'linewidth',2,'color',[0 0 1]);
for ii = 1:NH
    plot(ratios,gg_glm2(:,ii),'linewidth',1,'color',cols(ii,:));
end
xlabel('G_na/G_k');
ylabel('gain scaling:','FontSize',fontSizeLabel);
title('GLM fit to ','fontSize',fontSizeTitle)
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off


H = 16;
Na = 8;
K  = 8;

bps = squeeze(bitsPerSpike_test(:,Na,K,T,:,1))';

bps = nan(16,4);
for ii = 1:16
    for jj = 1:4
        bp = pseudo_R2_test(jj,:,:,T,ii,1);
        bps(ii,jj) = squeeze(nanmedian(bp(:)))';
    end
end

subplot(NR,NC,1);
for ii = 1:4
    plot(spkHistLengths,bps(:,ii),'linewidth',1,'color',sta_cs(ii,:));
    hold on
end
ylim([0 0.75]);
title('mean GLM predictive performance','fontSize',fontSizeTitle)
xlabel('spk history length (ms)');
ylabel('pesudo R^2','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off



bps = squeeze(bitsPerSpike_test(:,Na,K,:,H,1))';
bps = squeeze(nanmean(nanmean(bitsPerSpike_test(:,:,:,:,H,1))))';

bps = nan(5,4);
for ii = 1:5
    for jj = 1:4
        bp = pseudo_R2_test(jj,:,:,ii,H,1);
        bps(ii,jj) = squeeze(nanmean(bp(:)))';
    end
end
bps(bps < -1) = nan;
subplot(NR,NC,1+NC);
hold on
for ii = 1:4
    plot([1:4 5.5],bps(:,ii),'o-','linewidth',1,'color',sta_cs(ii,:));
end
ylim([-1 0.75]);
xlabel('training condition');
ylabel('pseudo R^2','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off







set(gcf,'PaperSize',[NC NR]*PlotSize,'PaperPosition',[0 0 NC NR]*PlotSize);

if(saveFigs)
    matfig2fyp(gcf,sprintf('%s/fig2_spkHistLength_raw.fyp', figDir));
end

%%
figure(3);
clf
H = 7;

NR = 2;
NC = 4;


tts_h = 1:size(spkHistBases{H},1);
tts_h = tts_h(tts_h > 20);

h_spks = spkHistBases{H}*reshape(glm_h_spk(1:size(spkHistBases{H},2),:,:,T,H),[],length(Gs)*length(Gs));
h_spks = h_spks(tts_h,order)';
h_spks = h_spks(~isnan(ratios),:);
h_spksn = h_spks./sqrt(sum(h_spks.^2,2));

k_stims = stimBasis*reshape(glm_k_stim(:,:,:,T,H),[],length(Gs)*length(Gs));
k_stims = k_stims(:,order)';
k_stims = k_stims(~isnan(ratios),:);
k_stimsn = k_stims./sqrt(sum(k_stims.^2,2));

[COEFF1, SCORE1, LATENT1, TSQUARED1, EXPLAINED1,MU1] = pca(h_spks,'NumComponents',3);
[COEFF2, SCORE2, LATENT2, TSQUARED2, EXPLAINED2,MU2] = pca(k_stims,'NumComponents',3);

PC1 = nan(length(Gs),length(Gs));
PC1(order(~isnan(ratios))) = SCORE1(:,1);
PC2 = nan(length(Gs),length(Gs));
PC2(order(~isnan(ratios))) = SCORE1(:,2);


subplot(NR,NC,1);
hold on
plot(COEFF2(:,1:2));
plot(MU2./sqrt(sum(MU2.^2)),'k');
ylabel('weight (normalized)','FontSize',fontSizeLabel);
xlabel('time (ms)','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off


subplot(NR,NC,3);
hold on

plot(ratios(~isnan(ratios)),SCORE2(:,1:2),'o');


xlabel('G_{Na}/G_{K}');
ylabel('PC weight','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off

subplot(NR,NC,4);
hold on
gg_hh_n = gg_hh(~isnan(ratios));
plot(gg_hh_n,SCORE2(:,1:2),'o');
xlabel('HH gain scaling: D_2');
ylabel('PC weight','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off

subplot(NR,NC,2);
hold on

m1 = nanmin(ratios);
m2 = nanmax(ratios);

ratios_n = ratios(~isnan(ratios));
plot(SCORE2(:,1),SCORE2(:,2),'ok');
% for ii = 1:size(SCORE2,1)
%     
%     cc = (ratios_n(ii)-m1)/(m2-m1)*0.9;
%     cc = [1 1 1]*cc;
%     
%     plot(SCORE2(ii,1),SCORE2(ii,2),'o','Color',cc);
% end
xlabel('PC 1','FontSize',fontSizeLabel);
ylabel('PC 2','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off


subplot(NR,NC,1*NC+1);
hold on
plot(tts_h,COEFF1(:,1:2));
plot(tts_h,MU1./sqrt(sum(MU1.^2)),'k');
ylabel('weight (normalized)','FontSize',fontSizeLabel);
xlabel('time (ms)','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off


subplot(NR,NC,1*NC+3);
hold on
plot(ratios(~isnan(ratios)),SCORE1(:,1:2),'o');
xlabel('G_{Na}/G_{K}');
ylabel('PC weight','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off

subplot(NR,NC,1*NC+4);
hold on
plot(gg_hh(~isnan(ratios)),SCORE1(:,1:2),'o');
xlabel('HH gain scaling: D_2');
ylabel('PC weight','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off

subplot(NR,NC,1*NC+2);
hold on
plot(SCORE1(:,1),SCORE1(:,2),'ok');
% for ii = 1:size(SCORE1,1)
%     
%     cc = (ratios_n(ii)-m1)/(m2-m1)*0.9;
%     cc = [1 1 1]*cc;
%     
%     plot(SCORE1(ii,1),SCORE1(ii,2),'o','Color',cc);
% end
xlabel('PC 1','FontSize',fontSizeLabel);
ylabel('PC 2','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off

set(gcf,'PaperSize',[NC NR]*PlotSize,'PaperPosition',[0 0 NC NR]*PlotSize);

if(saveFigs)
    matfig2fyp(gcf,sprintf('%s/fig2_pca_raw.fyp', figDir));
end
%%

figure(2);
clf
bps = squeeze(mean(pseudo_R2_test(5,:,:,T,H,1),1))';
subplot(NR,NC,0*NC+3);
hold on
dotPlot(Gs,Gs,bps,[],[0 1],[0 0 0]);
cb2 = colorbar();
ylabel(cb2,'pseudo R^2','FontSize',fontSizeLabel);
cb2.Visible = true;
xlim(lims);
ylim(lims);
xlabel('G_K');
ylabel('G_Na','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
set(gca,'XTick',600:400:2000,'YTick',600:400:2000);
axis square
hold off




bps = squeeze(pseudo_R2_test(4,:,:,2,H,1))';
bps = squeeze(mean(pseudo_R2_test(5,:,:,T2,H,1),1))';
subplot(NR,NC,1*NC+3);
hold on
dotPlot(Gs,Gs,bps,[],[0 1],[0 0 0]);
cb2 = colorbar();
ylabel(cb2,'pseudo R^2','FontSize',fontSizeLabel);
cb2.Visible = true;
xlim(lims);
ylim(lims);
xlabel('G_K');
ylabel('G_Na','FontSize',fontSizeLabel);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off
set(gcf,'PaperSize',[NC NR]*PlotSize,'PaperPosition',[0 0 NC NR]*PlotSize);
axis square
set(gca,'XTick',600:400:2000,'YTick',600:400:2000);
if(saveFigs)
    %matfig2fyp(gcf,sprintf('%s/fig1_GSsummary_raw.fyp', figDir));
    saveas(gcf,sprintf('%s/fig2_pseudoR2s_raw.pdf', figDir),'pdf');
end