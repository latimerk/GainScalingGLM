load('Results/Lundstrom/Lundstrom_phaseInfo.mat');
load('Results/Lundstrom/Lundstrom_timeInfo.mat');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromSims.mat', 'StimPeriods','StimLevels');
load('/home/latimerk/gitCode/GainScaling/Results/Lundstrom/LundstromGLMs_meta.mat', 'spkHistLengths');

spkHistLengths(1) = 10;


saveFigs = true;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 3;

hhs = [1 4 8 12];

xx = (1:30)./30-0.5/30;

legendObjs = cell(length(StimPeriods),1);
for ii = 1:length(StimPeriods)
    legendObjs{ii} = sprintf('T = %dms',StimPeriods(ii));
end

%%
figure(1);
clf;
for cc = 1:3
    subplot(2,3,cc);
    hold on
    
    plot(xx,fits_phase.cycleAvg_3AHP(:,:,cc+1));
    title(sprintf('\\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
    
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
    subplot(2,3,cc+3);
    hold on
    
    plot(xx,fits_phase.cycleAvg_Gain(:,:,cc+1));
    %title(sprintf('\\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
    
    if(cc == 1)
        ylabel('gain','FontSize',fontSizeLabel);
        xlabel('time (1/T)','FontSize',fontSizeLabel);
    end
    
    ylim([1 2]);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off
end

drawnow;
set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*2],'PaperPosition',[0 0 PlotSize*3 PlotSize*2]);
if(saveFigs)
    saveas(gcf,'Figs/FractionalDiff/HH_sine_avg.pdf');
end

figure(2);
clf;
for cc = 1:3
    subplot(2,3,cc);
    hold on
    
    plot(xx,fits_time.cycleAvg_3AHP(:,:,cc+1));
    title(sprintf('\\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
    
    if(cc == 1)
        ylabel('spks/s','FontSize',fontSizeLabel);
        xlabel('time (1/T)','FontSize',fontSizeLabel);
        legend(legendObjs,'FontSize',fontSizeLegend)
    end
    
    ylim([0 35]);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off
end

for cc = 1:3
    subplot(2,3,cc+3);
    hold on
    
    plot(xx,fits_time.cycleAvg_Gain(:,:,cc+1));
    %title(sprintf('\\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);
    
    if(cc == 1)
        ylabel('gain','FontSize',fontSizeLabel);
        xlabel('time (1/T)','FontSize',fontSizeLabel);
    end
    
    ylim([1 2]);
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    hold off
end

drawnow;
set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*2],'PaperPosition',[0 0 PlotSize*3 PlotSize*2]);
if(saveFigs)
    saveas(gcf,'Figs/FractionalDiff/HH_sq_avg.pdf');
end

%%

for tt = 1:3
    switch tt
        case 1
            ca = fits_phase.cycleAvg_glm_3AHP;
            ca2 = fits_time.cycleAvg_glm_sineTrain_3AHP;
            fName1 = 'GLM_sineTrain_sine_avg';
            fName2 = 'GLM_sineTrain_sq_avg';
        case 2
            ca = fits_phase.cycleAvg_glm_sqTrain_3AHP;
            ca2 = fits_time.cycleAvg_glm_3AHP;
            fName1 = 'GLM_squareTrain_sine_avg';
            fName2 = 'GLM_squareTrain_sq_avg';
        case 3
            ca = fits_phase.cycleAvg_glm_flatTrain_3AHP;
            ca2 = fits_time.cycleAvg_glm_flatTrain_3AHP;
            fName1 = 'GLM_flatTrain_sine_avg';
            fName2 = 'GLM_flatTrain_sq_avg';
    end
    
    figure(tt*10 + 1);
    clf;
    
    for hh = 1:length(hhs)
        hh_test = hhs(hh);
        for cc = 1:3
            subplot(length(hhs)+1,3,cc+(hh-1)*3);
            hold on

            plot(xx,ca(:,:,cc+1,hh_test));

            if(hh == 1 && cc == 1)
                ylabel('spks/s','FontSize',fontSizeLabel);
                xlabel('time (1/T)','FontSize',fontSizeLabel);
                
                title(sprintf('\\sigma = %.1f\nspk hist %dms',StimLevels(cc+1),spkHistLengths(hh_test)),'FontSize',fontSizeTitle);

                legend(legendObjs,'FontSize',fontSizeLegend)
            elseif(cc == 1)
                title(sprintf('spk hist %dms',spkHistLengths(hh_test)),'FontSize',fontSizeTitle);
            end

            ylim([3 24]);
            set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
            hold off
        end
    end
    for cc = 1:3
        subplot(length(hhs)+1,3,cc+length(hhs)*3);
        hold on

        plot(xx,fits_phase.cycleAvg_Gain(:,:,cc+1));
        %title(sprintf('\\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);

        if(cc == 1)
            ylabel('gain','FontSize',fontSizeLabel);
            xlabel('time (1/T)','FontSize',fontSizeLabel);
        end

        ylim([1 2]);
        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
        hold off
    end

    drawnow;
    set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*(length(hhs)+1)],'PaperPosition',[0 0 PlotSize*3 PlotSize*(length(hhs)+1)]);
    if(saveFigs)
        saveas(gcf,sprintf('Figs/FractionalDiff/%s.pdf',fName1));
    end

    figure(tt*10 + 2);
    clf;
    
    for hh = 1:length(hhs)
        hh_test = hhs(hh);
        for cc = 1:3
            subplot(length(hhs)+1,3,cc+(hh-1)*3);
            hold on

            plot(xx,ca2(:,:,cc+1,hh_test));
            title(sprintf('\\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);

                if(hh == 1 && cc == 1)
                    ylabel('spks/s','FontSize',fontSizeLabel);
                    xlabel('time (1/T)','FontSize',fontSizeLabel);

                    title(sprintf('\\sigma = %.1f\nspk hist %dms',StimLevels(cc+1),spkHistLengths(hh_test)),'FontSize',fontSizeTitle);

                    legend(legendObjs,'FontSize',fontSizeLegend)
                elseif(cc == 1)
                    title(sprintf('spk hist %dms',spkHistLengths(hh_test)),'FontSize',fontSizeTitle);
                end

            ylim([0 38]);
            set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
            hold off
        end
    end

    for cc = 1:3
        subplot(length(hhs)+1,3,cc+length(hhs)*3);
        hold on

        plot(xx,fits_time.cycleAvg_Gain(:,:,cc+1));
        %title(sprintf('\\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);

        if(cc == 1)
            ylabel('gain','FontSize',fontSizeLabel);
            xlabel('time (1/T)','FontSize',fontSizeLabel);
        end

        ylim([1 2]);
        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
        hold off
    end

    drawnow;
    set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*(length(hhs)+1)],'PaperPosition',[0 0 PlotSize*3 PlotSize*(length(hhs)+1)]);
    if(saveFigs)
        saveas(gcf,sprintf('Figs/FractionalDiff/%s.pdf',fName2));
    end

end
    
