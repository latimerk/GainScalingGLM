load Results/Mease/Mease_gainInfo.mat
load Results/Mease/MeaseGLMs.mat
load('/home/latimerk/gitCode/GainScaling/Results/Mease/gainDistance.mat')
figDir = 'Figs/GainScaling/';

addpath ~/FigureComposer/
addpath ~/FigureComposer/matlab/

saveFigs = false;


fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 2;
            
%%
Gs = 600:100:2000;

example_idx = [11 5;
                8 8;
                5 11];

sta_cs = [0.75 0.75 0.75;
          0.5  0.5  0.5;
          0.25 0.25 0.25;
          0    0    0];
      
sta_cs = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560];
      
sta_cs = [0    0    0;
          80   0  5;
          160  0  10;
          240 0 15]./255;
            
figure(1);
clf;

H = 16;
T = 5;

NE = size(example_idx,1);

NC = 7;

lims_p = [-2 4];
bins_p = fits_gain.bins >= lims_p(1) & fits_gain.bins <= lims_p(2);

for ii = 1:NE
    Na = example_idx(ii,1);
    K  = example_idx(ii,2);
    
    k_stim = stimBasis*glm_k_stim(:,Na,K,T,H);
    h_spk  = spkHistBases{H}*glm_h_spk(1:size(spkHistBases{H},2),Na,K,T,H);
    
    sta_glm   = squeeze(fits_gain.glm.sta_n(:,Na,K,:,T,H));
    p_spk_glm = squeeze(fits_gain.glm.p_spk(:,Na,K,:,T,H));
    p_all_glm = squeeze(fits_gain.glm.s_all(:,Na,K,:,T,H));
    p_all_glm = p_all_glm./sum(p_all_glm);
    f_spk_glm = p_spk_glm./p_all_glm;
    
    sta_hh = squeeze(fits_gain.hh.sta_n(:,Na,K,:));
    p_spk_hh = squeeze(fits_gain.hh.p_spk(:,Na,K,:));
    p_all_hh = squeeze(fits_gain.hh.s_all(:,Na,K,:));
    p_all_hh = p_all_hh./sum(p_all_hh);
    f_spk_hh = p_spk_hh./p_all_hh;
    
    subplot(NE,NC,1 + (ii-1)*NC);
    hold on
    plot([1 150],[0 0],'k:');
    plot(k_stim,'k','LineWidth',2);
    set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
    hold off
    
    subplot(NE,NC,2 + (ii-1)*NC);
    hold on
    plot([1 150],[0 0],'k:');
    for jj = 1:4
        plot(sta_glm(:,jj),'LineWidth',1,'Color',sta_cs(jj,:));
    end
    set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
    hold off
    
    subplot(NE,NC,3 + (ii-1)*NC);
    hold on
    for jj = 1:4
        plot(fits_gain.bins(bins_p), p_spk_glm(bins_p,jj),'LineWidth',1,'Color',sta_cs(jj,:));
    end
    xlim(lims_p);
    set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
    hold off
    
    subplot(NE,NC,4 + (ii-1)*NC);
    for jj = 1:4
        semilogy(fits_gain.bins(bins_p),    f_spk_glm(bins_p,jj),'LineWidth',1,'Color',sta_cs(jj,:));
        hold on
    end
    xlim(lims_p);
    set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
    hold off
    
    subplot(NE,NC,5 + (ii-1)*NC);
    hold on
    plot([1 150],[0 0],'k:');
    for jj = 1:4
        plot(sta_hh(:,jj),'LineWidth',1,'Color',sta_cs(jj,:));
    end
    set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
    hold off
    
    subplot(NE,NC,6 + (ii-1)*NC);
    hold on
    for jj = 1:4
        plot(fits_gain.bins(bins_p),    p_spk_hh(bins_p,jj),'LineWidth',1,'Color',sta_cs(jj,:));
    end
    xlim(lims_p);
    set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
    hold off
    
    
    subplot(NE,NC,7 + (ii-1)*NC);
    for jj = 1:4
        semilogy(fits_gain.bins(bins_p),    f_spk_hh(bins_p,jj),'LineWidth',1,'Color',sta_cs(jj,:));
        hold on
    end
    xlim(lims_p);
    set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
    hold off
end

%%
set(gcf,'PaperUnits','inches','PaperSize',[NC NE]*PlotSize,'PaperPosition',[0 0 NC NE]*PlotSize);
if(saveFigs)
    matfig2fyp(gcf,sprintf('%s/fig1_glmExample1_raw.fyp', figDir));
end