load('Results/Mease/HHexample.mat');
figDir = 'Figs/GainScaling/';

saveFigs = false;
tts_idx = (100*downsampleRate)+(1:(4000*downsampleRate));
tts = (1:length(tts_idx))./downsampleRate;

xlims = [0 tts(end)];

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 3;


figure(1);
clf

for ii = 1:4
    subplot(2,4,ii);
    
    
    plot(tts,y(tts_idx)*sigs(ii) + stim_dc,'k');
    
    title(sprintf('\\sigma = %.1f',StimLevels(ii)),'FontSize',fontSizeTitle);
    if(ii == 1)
        ylabel('\mu A/cm^2','FontSize',fontSizeLabel);
        xlabel('time (ms)','FontSize',fontSizeLabel);
    end
    ylim([-6 6]);
    xlim(xlims);
    
    set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
    
    subplot(2,4,ii+4)
    plot(tts,squeeze(V(1,ii,tts_idx)),'k');
    
    if(ii == 1)
        ylabel('mV','FontSize',fontSizeLabel);
        xlabel('time (ms)','FontSize',fontSizeLabel);
    end
    ylim([-90 40]);
    xlim(xlims);
    
    set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
end

set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*4 PlotSize*2],'PaperPosition',[0 0 PlotSize*4 PlotSize*2]);
if(saveFigs)
    saveas(gcf,sprintf('%s/exampleHH.pdf', figDir));
end