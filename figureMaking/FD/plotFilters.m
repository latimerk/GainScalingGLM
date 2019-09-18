load('Results/Lundstrom/LundstromGLMs_meta.mat');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat', 'spkHistLengths');

spkHistLengths(1) = 10;


saveFigs = false;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 3;

hhs = [1 6 8 10 12];
%%
cc = 4;
for hh = 1:length(hhs)
    hh_test = hhs(hh);
    figure(hh);
    clf
    subplot(2,2,1)
    hold on
    plot([1 size(stimBasis,1)],[0 0]','k--');
    plot(stimBasis*squeeze(glm_k_stim(:,end,cc,hh_test,5)),'Color',[0    0.4470    0.7410]);
    plot(stimBasis*squeeze(glm_k_stim(:,end,cc,hh_test,2)),'Color',[0.8500    0.3250    0.0980]);
    plot(stimBasis*squeeze(glm_k_stim(:,end,1  ,hh_test,5)),'Color',[0.6 0.6 0.6]);
    xlabel('time (ms)','FontSize',fontSizeLabel);
    ylabel('filter (a.u.)','FontSize',fontSizeLabel);
    
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    title(sprintf('stimulus filter\n\\sigma = %.1f',StimLevels(cc)),'fontSize',fontSizeTitle);
    xlim([1 size(stimBasis,1)]);
    
    legend({'GLM_{sine}','GLM_{square}','GLM_{flat}'},'FontSize',fontSizeLegend)
    
    hold off
    
    subplot(2,2,2)
    hold on
    
    plot([1 size(spkHistBases{hh_test},1)],[0 0]','k--');
    plot(spkHistBases{hh_test}*squeeze(glm_h_spk(1:size(spkHistBases{hh_test},2),end,cc,hh_test,5)),'Color',[0    0.4470    0.7410]);
    plot(spkHistBases{hh_test}*squeeze(glm_h_spk(1:size(spkHistBases{hh_test},2),end,cc,hh_test,2)),'Color',[0.8500    0.3250    0.0980]);
    plot(spkHistBases{hh_test}*squeeze(glm_h_spk(1:size(spkHistBases{hh_test},2),end,1  ,hh_test,5)),'Color',[0.6 0.6 0.6]);
    
    xlabel('time (ms)','FontSize',fontSizeLabel);
    ylabel('filter (a.u.)','FontSize',fontSizeLabel);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    title(sprintf('spike history filter (%d ms)',spkHistLengths(hh_test)),'fontSize',fontSizeTitle);
    xlim([1 size(spkHistBases{hh_test},1)]);
    hold off
    
    
    subplot(2,2,3)
    hold on
    plot([1 size(stimBasis,1)],[1 1]','k--');
    plot(exp(stimBasis*squeeze(glm_k_stim(:,end,cc,hh_test,5))),'Color',[0    0.4470    0.7410]);
    plot(exp(stimBasis*squeeze(glm_k_stim(:,end,cc,hh_test,2))),'Color',[0.8500    0.3250    0.0980]);
    plot(exp(stimBasis*squeeze(glm_k_stim(:,end,1  ,hh_test,5))),'Color',[0.6 0.6 0.6]);
    xlabel('time (ms)','FontSize',fontSizeLabel);
    ylabel('gain','FontSize',fontSizeLabel);
    
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    %title('stimulus filter','fontSize',fontSizeTitle);
    xlim([1 size(stimBasis,1)]);
    hold off
    
    subplot(2,2,4)
    
    semilogx([1 size(spkHistBases{hh_test},1)],[1 1]','k--');
    hold on
    plot(exp(spkHistBases{hh_test}*squeeze(glm_h_spk(1:size(spkHistBases{hh_test},2),end,cc,hh_test,5))),'Color',[0    0.4470    0.7410]);
    plot(exp(spkHistBases{hh_test}*squeeze(glm_h_spk(1:size(spkHistBases{hh_test},2),end,cc,hh_test,2))),'Color',[0.8500    0.3250    0.0980]);
    plot(exp(spkHistBases{hh_test}*squeeze(glm_h_spk(1:size(spkHistBases{hh_test},2),end,1  ,hh_test,5))),'Color',[0.6 0.6 0.6]);
    xlim([1 size(spkHistBases{hh_test},1)]);
    xlabel('time (ms)','FontSize',fontSizeLabel);
    ylabel('gain','FontSize',fontSizeLabel);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    %title(sprintf('spike history filter (%d ms)',spkHistLengths(hh_test)),'fontSize',fontSizeTitle);
    hold off
    
    
    set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*2 PlotSize*2],'PaperPosition',[0 0 PlotSize*2 PlotSize*2]);
    if(saveFigs)
        saveas(gcf,sprintf('Figs/FractionalDiff/glm_filters_%d_sig%d.pdf',hh_test,cc));
    end
    
end

%%
cc = 4;
figure(10);
clf
subplot(1,3,1)
hold on

ccs = linspace(0.1,1,length(hhs));
for hh = 1:length(hhs)
    hh_test = hhs(hh);
    plot(stimBasis*squeeze(glm_k_stim(:,end,cc,hh_test,5)),'Color',ccs(hh)*[0    0.4470    0.7410]);
end
    

title(sprintf('stimulus filters (\\sigma = %.1f)\nGLM_{sine}',StimLevels(cc)),'fontSize',fontSizeTitle);

xlabel('time (ms)','FontSize',fontSizeLabel);
ylabel('filter (a.u.)','FontSize',fontSizeLabel);

set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
%title('stimulus filter','fontSize',fontSizeTitle);
xlim([1 size(stimBasis,1)]);

legendObjs = cell(length(hhs),1);
for hh = 1:length(hhs)
    legendObjs{hh} = sprintf('spk hist %dms',spkHistLengths(hhs(hh)));
end
legend(legendObjs);

hold off



subplot(1,3,2)
hold on

for hh = 1:length(hhs)
    hh_test = hhs(hh);
    plot(stimBasis*squeeze(glm_k_stim(:,end,cc,hh_test,2)),'Color',ccs(hh)*[0.8500    0.3250    0.0980]);
end

xlabel('time (ms)','FontSize',fontSizeLabel);
ylabel('filter (a.u.)','FontSize',fontSizeLabel);
title('GLM_{square}','fontSize',fontSizeTitle);

legendObjs = cell(length(hhs),1);
for hh = 1:length(hhs)
    legendObjs{hh} = sprintf('spk hist %dms',spkHistLengths(hhs(hh)));
end
legend(legendObjs);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
%title('stimulus filter','fontSize',fontSizeTitle);
xlim([1 size(stimBasis,1)]);
hold off


subplot(1,3,3)
hold on

for hh = 1:length(hhs)
    hh_test = hhs(hh);
    plot(stimBasis*squeeze(glm_k_stim(:,end,1,hh_test,2)),'Color',ccs(hh)*[0.6 0.6 0.6]);
end

xlabel('time (ms)','FontSize',fontSizeLabel);
ylabel('filter (a.u.)','FontSize',fontSizeLabel);
title('GLM_{flat}','fontSize',fontSizeTitle);

legendObjs = cell(length(hhs),1);
for hh = 1:length(hhs)
    legendObjs{hh} = sprintf('spk hist %dms',spkHistLengths(hhs(hh)));
end
legend(legendObjs);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
%title('stimulus filter','fontSize',fontSizeTitle);
xlim([1 size(stimBasis,1)]);
hold off

set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*1],'PaperPosition',[0 0 PlotSize*3 PlotSize*1]);
if(saveFigs)
    saveas(gcf,sprintf('Figs/FractionalDiff/glm_stimFilters_sig%d.pdf',cc));
end
