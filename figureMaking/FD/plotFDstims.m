load('Results/Lundstrom/HHexample.mat');
figDir = 'Figs/FractionalDiff/';

saveFigs = false;
tts_idx = (2000*downsampleRate)+(1:(4000*downsampleRate));
tts = (1:length(tts_idx))./downsampleRate;

xlims = [0 tts(end)];

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 3;


figure(1);
clf

subplot(2,2,1);


plot(tts,y_sine(tts_idx),'k');

title(sprintf('sine modulated, \\sigma = %.1f',StimLevel),'FontSize',fontSizeTitle);

ylabel('{\mu}A/cm^2','FontSize',fontSizeLabel);
xlabel('time (ms)','FontSize',fontSizeLabel);

ylim([-25 25]);
xlim(xlims);

set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

subplot(2,2,3)
plot(tts,squeeze(V_sine(1,1,tts_idx)),'k');

ylabel('mV','FontSize',fontSizeLabel);
xlabel('time (ms)','FontSize',fontSizeLabel);

ylim([-95 50]);
xlim(xlims);

set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);




subplot(2,2,2);


plot(tts,y_sq(tts_idx),'k');

title(sprintf('square modulated, \\sigma = %.1f',StimLevel),'FontSize',fontSizeTitle);


ylim([-25 25]);
xlim(xlims);

set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

subplot(2,2,4)
plot(tts,squeeze(V_sq(1,1,tts_idx)),'k');


ylim([-95 50]);
xlim(xlims);

set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);


set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*2],'PaperPosition',[0 0 PlotSize*3 PlotSize*2]);
if(saveFigs)
    saveas(gcf,sprintf('%s/exampleHH.pdf', figDir));
end