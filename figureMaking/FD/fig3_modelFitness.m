
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/glmFitMetrics_4.mat')
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
figure(3);
clf;
dd = 1e3;

NR = 4;
NC = 3;

Hs = 1:26;%[1 5 8 15 26];
NH = length(Hs);
cols = repmat(linspace(0.9,0,NH)',[1 3]);
cols2 = repmat(linspace(0.9,0,26)',[1 3]);
ffs = freqs <= 32 & (1./freqs <= 64);
ffs = (1./freqs <= 64) & (1./freqs >= 0.01);

L = 4;
L2 = L;

colorSineFit = [0 0 0];
colorSquareFit = [0 0 0];%[0    0.4470    0.7410];
colorFlatFit = [    0.8500    0.3250    0.0980];

lims= [0.5 0.65];

% pseudo r2: sine noise
subplot(NR,NC,1);
semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,1,:,1)),'-','Color',colorSineFit);
hold on
% semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,1,:,2)),'-','Color',colorSquareFit);

semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,2,:,1)),'--','Color',colorSineFit);
% semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,2,:,2)),'--','Color',colorSquareFit);

% semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,1,:,3)),'-','Color',colorFlatFit);
xlim([10 16e3]./dd);
if(dd == 1e3)
    xlabel('spike history length (s)','FontSize',fontSizeLabel);
else
    xlabel('spike history length (ms)','FontSize',fontSizeLabel);
end
ylim(lims);
ylabel('pseudo-','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off


%pseudo r2: square noise
subplot(NR,NC,2);
semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,1,:,2)),'-','Color',colorSquareFit);
hold on
semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,2,:,2)),'--','Color',colorSquareFit);
% semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,2,:,3)),'-','Color',colorFlatFit);
xlim([10 16e3]./dd);
if(dd == 1e3)
    xlabel('spike history length (s)','FontSize',fontSizeLabel);
else
    xlabel('spike history length (ms)','FontSize',fontSizeLabel);
end
ylim(lims);
ylabel('pseudo-','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off

% pseudo r2: flat train
subplot(NR,NC,3);
semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,1,:,3)),'-','Color',colorFlatFit);
hold on
semilogx(spkHistLengths./dd,squeeze(pseudo_R2_test(L,2,:,3)),'--','Color',colorFlatFit);
xlim([10 16e3]./dd);
if(dd == 1e3)
    xlabel('spike history length (s)','FontSize',fontSizeLabel);
else
    xlabel('spike history length (ms)','FontSize',fontSizeLabel);
end
ylim([-20 0.65]);
ylabel('pseudo-','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off

% coherence: sine fit
subplot(NR,NC,4);
for hh = 1:NH
    semilogx(1./freqs(ffs),squeeze(mean(c_test(ffs,L2,1,Hs(hh),1),2)),'Color',cols(hh,:))
    hold on
end
% xlim([min(freqs(ffs)) max(freqs(ffs))]);
ylabel('coherence','FontSize',fontSizeLabel);
xlabel('period (s)','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off


% coherence: square fit
subplot(NR,NC,5);
for hh = 1:NH
    semilogx(1./freqs(ffs),squeeze(mean(c_test(ffs,L2,1,Hs(hh),2),2)),'Color',cols(hh,:))
    hold on
end
% xlim([min(freqs(ffs)) max(freqs(ffs))]);
ylabel('coherence','FontSize',fontSizeLabel);
xlabel('period (s)','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off

% % coherence: flat fit
subplot(NR,NC,6);
for hh = 1:NH
    semilogx(1./freqs(ffs),squeeze(mean(c_test(ffs,L2,1,Hs(hh),3),2)),'Color',cols(hh,:))
    hold on
end
% xlim([min(freqs(ffs)) max(freqs(ffs))]);
ylabel('coherence','FontSize',fontSizeLabel);
xlabel('period (s)','FontSize',fontSizeLabel);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
hold off


for CC = 1:3
    %%
    hs = zeros(size(spkHistBases{16},1),16);
    ks = zeros(size(stimBasis,1),16);

    for hh = 1:26

        if(CC == 1)
            CC_c = 5;
            Lc = L;
        elseif(CC == 2)
            CC_c = 2;
            Lc = L;
        else
            CC_c = 5;
            Lc = 1;
        end
        hs(1:size(spkHistBases{hh} ,1),hh) = spkHistBases{hh }*glm_h_spk(1:size(spkHistBases{hh },2),end,Lc,hh ,CC_c);
        ks(:,hh) = stimBasis*glm_k_stim(:,end,Lc,hh,CC_c);
    end


    subplot(NR,NC,6+CC);
    for hh = 1:26
        semilogx((1:size(hs,1))./dd,cumsum(hs(:,hh)),'Color',cols2(hh,:));
        hold on
    end
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off


    subplot(NR,NC,6+NC+CC);
    for hh = 1:26
        plot((ks(:,hh)),'Color',cols2(hh,:));
        hold on
    end
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off
end

set(gcf,'PaperUnits','inches','PaperSize',[NC NR]*PlotSize,'PaperPosition',[0 0 NC NR]*PlotSize);
if(saveFigs)
    matfig2fyp(gcf,sprintf('%s/fig3_modelFitness_raw.fyp', figDir));
end
