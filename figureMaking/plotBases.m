

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;

addpath GLMtools/
[~,stimBasis] = getStimBasis();
[~,sb_GS] = getSpkHistBasis((150-1)/1e3,15);
[~,sb_FD] = getSpkHistBasis((16e3-1)/1e3,25);


figure(1);
clf

subplot(3,1,1)
plot(stimBasis);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);     

subplot(3,1,2)
plot(sb_GS);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);     

subplot(3,1,3)
plot(sb_FD);
set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);     


set(gcf,'PaperSize',[3 3],'PaperPosition',[0 0 3 3]);

addpath ~/FigureComposer/
addpath ~/FigureComposer/matlab/