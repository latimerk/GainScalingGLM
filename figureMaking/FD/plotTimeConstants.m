
load('Results/Lundstrom/Lundstrom_timeInfo.mat');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat', 'spkHistLengths');
spkHistLengths(1) = 10;
saveFigs = false;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 3;

hhs = [6 8 12];

figure(1);
clf;

for hh = 1:length(hhs)
    for ii = 2:4
        subplot(length(hhs),3,ii-1 + 3*(hh-1))
        hh_test = hhs(hh);
        hold on

        
        
        
        %p5 = plot(StimPeriods,(fits_time.tau_Default(:,ii,1)),'k^:','LineWidth',1.5);
        p3 = plot(StimPeriods,(fits_time.tau_glm_sineTrain_3AHP(:,ii,hh_test,1)),'^-','LineWidth',1.5,'Color',[0    0.4470    0.7410]);
        p2 = plot(StimPeriods,(fits_time.tau_glm_3AHP(:,ii,hh_test,1)),'^-','LineWidth',1.5,'Color',[ 0.8500    0.3250    0.0980]);
        p4 = plot(StimPeriods,(fits_time.tau_glm_flatTrain_3AHP(:,ii,hh_test,1)),'^-','Color',[0.6 0.6 0.6],'LineWidth',1.5);
        p1 = plot(StimPeriods,(fits_time.tau_3AHP(:,ii,1)),'k^-','LineWidth',2);

        %plot(StimPeriods,(fits_time.tau_Default(:,ii,2)),'v:','Color',p5.Color,'LineWidth',p5.LineWidth);
        plot(StimPeriods,(fits_time.tau_glm_sineTrain_3AHP(:,ii,hh_test,2)),'v-','Color',p3.Color,'LineWidth',p3.LineWidth);
        plot(StimPeriods,(fits_time.tau_glm_3AHP(:,ii,hh_test,2)),'v-','Color',p2.Color,'LineWidth',p2.LineWidth);
        plot(StimPeriods,(fits_time.tau_glm_flatTrain_3AHP(:,ii,hh_test,2)),'v-','Color',p4.Color,'LineWidth',p4.LineWidth);
        plot(StimPeriods,(fits_time.tau_3AHP(:,ii,2)),'v-','Color',p1.Color,'LineWidth',p1.LineWidth);
        
        ylabel('Time constant \tau (s)','FontSize',fontSizeLabel);
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
        ylim([-0.5 4]);

        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
        hold off
    end
end

drawnow;
set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*length(hhs)],'PaperPosition',[0 0 PlotSize*3 PlotSize*length(hhs)]);
if(saveFigs)
    saveas(gcf,'Figs/FractionalDiff/boxTimeConstants.pdf');
end