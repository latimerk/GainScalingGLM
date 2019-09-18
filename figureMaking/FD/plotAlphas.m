load('Results/Lundstrom/Lundstrom_phaseInfo.mat');
load('Results/Lundstrom/Lundstrom_timeInfo.mat');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat', 'spkHistLengths');
load('Results/Lundstrom/alphas.mat');
addpath Utils/
spkHistLengths(1) = 10;

HH = length(spkHistLengths);
saveFigs = false;
NB = 30;


%%

hhs = 1:HH;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 3;

figure(1);
clf;
for xx = 1:3
    switch xx
        case 1
            bb = alphas.filter;
        case 2
            bb = alphas.gain;
        case 3
            bb = alphas.phase;
    end
    
    for cc = 1:3
        subplot(3,3,cc + (xx-1)*3);
        
        semilogx(spkHistLengths(hhs),ones(size(spkHistLengths(hhs)))*bb.hh_def(cc),':','LineWidth',2,'Color',[0 0 0]);
        hold on
        semilogx(spkHistLengths(hhs),bb.sine(hhs,cc)  ,'LineWidth',1.5,'Color',[0    0.4470    0.7410]);
        plot(spkHistLengths(hhs),bb.square(hhs,cc),'LineWidth',1.5,'Color',[0.8500    0.3250    0.0980]);
        plot(spkHistLengths(hhs),bb.flat(hhs,cc)  ,'LineWidth',1.5,'Color',[0.6 0.6 0.6]);
        plot(spkHistLengths(hhs),ones(size(spkHistLengths(hhs)))*bb.hh(cc),'LineWidth',2,'Color',[0 0 0]);
        
        if(xx == 1)
            if(cc == 1)
                title(sprintf('Box Noise: Filter\n\\sigma=%.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
            else
                title(sprintf('\\sigma=%.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
            end
        elseif(cc == 1)
            if(xx == 2)
                title('Sine Noise: Gain','FontSize',fontSizeTitle);
            elseif(xx == 3)
                title('Sine Noise: Phase','FontSize',fontSizeTitle);
            end
        end
        
        if(xx == 1 && cc == 1)
            cs = get(gca,'Children');
            legendObs = {'GLM_{sine}','GLM_{square}','GLM_{flat}','Hodgkin-Huxley plus AHP','Hodgkin-Huxley no AHP'};
            legend(cs([2 3 4 1 5]),legendObs,'Location','SouthEast','FontSize',fontSizeLegend);
        end
        
        
        ylim([-0.01 0.3]);
        xlim([0 max(spkHistLengths(hhs))]);
        ylabel('fractional order \alpha','FontSize',fontSizeLabel);
        xlabel('spk hist length (ms)','FontSize',fontSizeLabel);
        
        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
        hold off
    end
end

drawnow;
set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*3],'PaperPosition',[0 0 PlotSize*3 PlotSize*3]);
if(saveFigs)
    saveas(gcf,'Figs/FractionalDiff/estimatedAlphas.pdf');
end



