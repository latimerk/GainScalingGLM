load('Results/Lundstrom/Lundstrom_phaseInfo.mat');
load('Results/Lundstrom/Lundstrom_timeInfo.mat');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat', 'spkHistLengths');

spkHistLengths(1) = 10;

addpath ~/FigureComposer/
addpath ~/FigureComposer/matlab/


saveFigs = true;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 2;

figDir = 'Figs/FractionalDiff/';

%%

T = 4;
H = 26;

NR = 4;
NC = 3;

colors = [0       2*16+2   2*16+11;
          0       6*16+2   8*16 + 0;
          0       8*16+8   10*16+10;
          0       10*16+10 13*16+4;
          0       12*16+12 15*16+15;
          2*16+10 13*16+4  15*16+15;
          5*16+5  13*16+13 15*16+15]./255;
colors = [ 0 0 0;
           1*16+12    2*16+ 4    2*16+2;
           3*16+ 7    4*16+ 8    4*16+5; 
           5*16+ 3    6*16+12    6*16+7; 
           6*16+15    9*16+ 1    8*16+10; 
           9*16+ 3   10*16+12   10*16+7; 
          10*16+ 7   12*16+ 8   12*16+4]./255;
      
      
% colors = get(gca,'colororder');
% colors = repmat(linspace(0.8,0,7)',[1 3]);


      
NF = 7;
figure(1);
clf


xx = (1:30)./30-0.5/30;
legendObjs = cell(length(StimPeriods),1);
for ii = 1:length(StimPeriods)
    legendObjs{ii} = sprintf('T = %dms',StimPeriods(ii));
end


for cc = 1:3
    subplot(NR,NC,cc);
    hold on
    
    for ff = 1:NF
        plot(xx,fits_phase.cycleAvg_3AHP(:,ff,cc+1),'Color',colors(ff,:));
    end
    title(sprintf('Hodgkin-Huxley: %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
    
    if(cc == 1)
        ylabel('spks/s','FontSize',fontSizeLabel);
        xlabel('time (1/T)','FontSize',fontSizeLabel);
    
        legend(legendObjs,'FontSize',fontSizeLegend)
    end
    
    ylim([3 24]);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off
end



for cc = 1:3
    subplot(NR,NC,NC*1+cc);
    hold on
    
    for ff = 1:NF
        plot(xx,fits_phase.cycleAvg_glm_3AHP(:,ff,cc+1,H),'Color',colors(ff,:));
    end
    title(sprintf('GLM: %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
    
    if(cc == 1)
        ylabel('spks/s','FontSize',fontSizeLabel);
        xlabel('time (1/T)','FontSize',fontSizeLabel);
    end
    
    ylim([3 24]);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off
end



for cc = 1:3
    subplot(NR,NC,cc+NC*2);
    hold on
    
    for ff = 1:NF
        plot(xx,fits_time.cycleAvg_3AHP(:,ff,cc+1),'Color',colors(ff,:));
    end
    title(sprintf('Hodgkin-Huxley: %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
    
    if(cc == 1)
        ylabel('spks/s','FontSize',fontSizeLabel);
        xlabel('time (1/T)','FontSize',fontSizeLabel);
    end
    
    ylim([0 33]);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off
end



for cc = 1:3
    subplot(NR,NC,cc+NC*3);
    hold on
    
    for ff = 1:NF
        plot(xx,fits_time.cycleAvg_glm_3AHP(:,ff,cc+1,H),'Color',colors(ff,:));
    end
    title(sprintf('GLM: %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
    
    if(cc == 1)
        ylabel('spks/s','FontSize',fontSizeLabel);
        xlabel('time (1/T)','FontSize',fontSizeLabel);
    end
    
    ylim([0 33]);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off
end


set(gcf,'PaperUnits','inches','PaperSize',[NC NR]*PlotSize,'PaperPosition',[0 0 NC NR]*PlotSize);
if(saveFigs)
    matfig2fyp(gcf,sprintf('%s/fig1_cycleAvgs_raw.fyp', figDir));
end