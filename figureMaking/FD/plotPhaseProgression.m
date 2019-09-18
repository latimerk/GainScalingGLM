
load('Results/Lundstrom/Lundstrom_phaseInfo.mat');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat', 'spkHistLengths');
saveFigs = false;
spkHistLengths(1) = 10;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 3;

hhs = [1 4 8 12];

figure(1);
clf;
for hh = 1:length(hhs)
    for ii = 2:4
        subplot(length(hhs),3,ii-1 + 3*(hh-1))
        hh_test = hhs(hh);
        hold on
        plot(StimPeriods,rad2deg(fits_phase.phase_Default(:,ii)-fits_phase.phase_Gain(:,ii)),'ko:','LineWidth',2)
        plot(StimPeriods,rad2deg(fits_phase.phase_glm_3AHP(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'o-','LineWidth',1.5)
        plot(StimPeriods,rad2deg(fits_phase.phase_glm_sqTrain_3AHP(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'o-','LineWidth',1.5)
        plot(StimPeriods,rad2deg(fits_phase.phase_glm_flatTrain_3AHP(:,ii,hh_test)-fits_phase.phase_Gain(:,ii)),'o-','LineWidth',1.5,'Color',[0.6 0.6 0.6])
        plot(StimPeriods,rad2deg(fits_phase.phase_3AHP(:,ii)-fits_phase.phase_Gain(:,ii)),'ko-','LineWidth',2)
        
        
        
        ylabel('Phase lead \psi (deg)','FontSize',fontSizeLabel);
        xlabel('Period T (ms)','FontSize',fontSizeLabel);
        
        if(hh == 1)
            if(ii == 2)
                title(sprintf('\\sigma = %.1f\n%d ms spk hist',StimLevels(ii),   spkHistLengths(hh_test)),'FontSize',fontSizeTitle);
            else
                title(sprintf('\\sigma = %.1f',StimLevels(ii)),'FontSize',fontSizeTitle);
            end
        elseif(ii == 2)
            title(sprintf('%d ms spk hist',spkHistLengths(hh_test)),'FontSize',fontSizeTitle);
        end
        
        xlim([0 64]);
        ylim([-4 24]);

        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
        hold off
    end
end


drawnow;
set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*length(hhs)],'PaperPosition',[0 0 PlotSize*3 PlotSize*length(hhs)]);
if(saveFigs)
    saveas(gcf,'Figs/FractionalDiff/phaseLeads.pdf');
end