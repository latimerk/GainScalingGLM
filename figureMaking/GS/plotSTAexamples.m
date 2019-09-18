load Results/Mease/Mease_gainInfo.mat
load Results/Mease/MeaseGLMs.mat
load('/home/latimerk/gitCode/GainScaling/Results/Mease/MeaseSims.mat', 'stim_mu','xs','frFunctions', 'G_mat_info_paper')

saveFigs = true;

spkHistLengths(1) = 10;
fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 3;
HH = length(spkHistLengths);

Gs = 600:100:2000;

figDir = 'Figs/GainScaling/';

ggs = [5 1;
       5 5;
       1 10];

%%

for gg = 1:size(ggs,1)
    G_Na = ggs(gg,1);
    G_K = ggs(gg,2);

    figure(gg);
    clf;

    subplot(1+length(hhs),2,1)
    hold on

    plot(squeeze(fits_gain.hh.sta_n(:,G_Na,G_K,:)));

    title(sprintf('Hodgkin-Huxley\nG_{Na} = %d, G_{K} = %d',Gs(G_Na),Gs(G_K)),'FontSize',fontSizeTitle);
    xlabel('time (ms)','FontSize',fontSizeLabel);
    ylabel('normalized sta','FontSize',fontSizeLabel);

    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);

    subplot(1+length(hhs),2,2)
    hold on

    % plot(squeeze(fits_gain.hh.sta_n(:,G_Na,G_K,:)));
    plot(fits_gain.bins,   squeeze(fits_gain.hh.p_spk(:,G_Na,G_K,:)));

    % title(sprintf('Hodgkin-Huxley\nG_{Na} = %d, G_{K} = %d',Gs(G_Na),Gs(G_K)),'FontSize',fontSizeTitle);
    xlabel('stim. \hat{s}','FontSize',fontSizeLabel);
    ylabel('spike-triggered stim. probability','FontSize',fontSizeLabel);

    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);


    hhs = [1 3 5 10];
    for hh = 1:length(hhs)
        hh_test = hhs(hh);

        subplot(1+length(hhs),2,2*(hh)+1)
        hold on

        plot(squeeze(fits_gain.glm.sta_n(:,G_Na,G_K,:,end,hh_test)));

        title(sprintf('GLM\nspk hist %d ms',spkHistLengths(hh_test)),'FontSize',fontSizeTitle);
        xlabel('time (ms)','FontSize',fontSizeLabel);
        ylabel('normalized sta','FontSize',fontSizeLabel);

        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);

        subplot(1+length(hhs),2,2*(hh)+2)
        hold on

        % plot(squeeze(fits_gain.hh.sta_n(:,G_Na,G_K,:)));
        plot(fits_gain.bins,   squeeze(fits_gain.glm.p_spk(:,G_Na,G_K,:,end,hh_test)));

        % title(sprintf('Hodgkin-Huxley\nG_{Na} = %d, G_{K} = %d',Gs(G_Na),Gs(G_K)),'FontSize',fontSizeTitle);
        xlabel('stim. \hat{s}','FontSize',fontSizeLabel);
        ylabel('spike-triggered stim. probability','FontSize',fontSizeLabel);

        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    end

    set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*2 PlotSize*(length(hhs)+1)],'PaperPosition',[0 0 PlotSize*2 PlotSize*(length(hhs)+1)]);
    if(saveFigs)
        saveas(gcf,sprintf('%s/GLM_sta_p_spk_GNa%d_GK%d.pdf', figDir,Gs(G_Na),Gs(G_K)));
    end
end

%%
figure(2);
clf; hold on

hhs = 1:HH;

for gg = 1:size(ggs,1)
    G_Na = ggs(gg,1);
    G_K = ggs(gg,2);

    subplot(size(ggs,1),2,1+(gg-1)*2)
    hold on


    ccs = linspace(0.8,0,length(hhs));
    legendObjs = cell(length(hhs),1);
    for ii = 1:length(hhs)
        plot(stimBasis* squeeze(glm_k_stim(:,G_Na,G_K,end,hhs(ii))),'Color',[1 1 1]*ccs(ii));
        legendObjs{ii} = sprintf('%dms',spkHistLengths(hhs(ii)));
    end
    if(gg == 1)
        legend(legendObjs,'FontSize',fontSizeLegend);
    end
    xlabel('time (ms)','FontSize',fontSizeLabel);
    ylabel('stimulus filter','FontSize',fontSizeLabel);
    if(gg == 1)
        title(sprintf('stimulus filters\n over spike history lengths\n G_{Na} = %d, G_{K} = %d',Gs(G_Na),Gs(G_K)),'FontSize',fontSizeTitle)
    else
        title(sprintf('G_{Na} = %d, G_{K} = %d',Gs(G_Na),Gs(G_K)),'FontSize',fontSizeTitle)
    end
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    
    
    
    subplot(size(ggs,1),2,2+(gg-1)*2)
    hold on


    ccs = linspace(0.8,0,length(hhs));
    legendObjs = cell(length(hhs),1);
    for ii = 1:length(hhs)
        plot(exp(spkHistBases{ii}* squeeze(glm_h_spk(1:size(spkHistBases{ii},2),G_Na,G_K,end,hhs(ii)))),'Color',[1 1 1]*ccs(ii));
        legendObjs{ii} = sprintf('%dms',spkHistLengths(hhs(ii)));
    end

    xlabel('time (ms)','FontSize',fontSizeLabel);
    ylabel('spk hist gain','FontSize',fontSizeLabel);

    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);

end

set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*2 PlotSize*length(ggs)],'PaperPosition',[0 0 PlotSize*2 PlotSize*length(ggs)]);
if(saveFigs)
    saveas(gcf,sprintf('%s/GLM_filterExamples.pdf', figDir));
end