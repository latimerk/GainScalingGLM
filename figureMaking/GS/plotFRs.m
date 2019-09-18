load Results/Mease/Mease_gainInfo.mat
load Results/Mease/MeaseGLMs.mat


addpath Utils
addpath Utils/EMD/
addpath figureMaking/

saveFigs = true;


fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 3;

Gs = 500 + (1:15)*100;

xticks = 600:400:2000;
yticks = 600:400:2000;

xrange = [Gs(1)-50 Gs(end)+50];
yrange = [Gs(1)-50 Gs(end)+50];

tts = [1 5];
hhs = [1:5 10];

hh_fr_range = [0 ceil(nanmax(fits_gain.hh.r_bar(:)))];

glm_fr_range = [0 ceil(nanmax(fits_gain.glm.r_bar(:)))];

all_fr_range = [0 max(hh_fr_range(2),glm_fr_range(2))];

TrainStrings = cell(length(StimLevels)+1,1);
for ii = 1:length(StimLevels)
    TrainStrings{ii} = sprintf('\\sigma = %.1f',StimLevels(ii));
end
TrainStrings{end} = 'all \sigma';

figDir = 'Figs/GainScaling/';
if(~isfolder(figDir))
    mkdir(figDir);
end
%%

figure(1);
clf
for cc = 1:4
    subplot(1,4,cc);
    dotPlot(Gs,Gs,fits_gain.hh.r_bar(:,:,cc)',[],all_fr_range);

    if(cc == 1)
        xlabel('G_K [pS/{\mu}m^2]','FontSize',fontSizeLabel);
        ylabel('G_{Na} [pS/{\mu}m^2]','FontSize',fontSizeLabel);
    end
    title(sprintf('stimulus: \\sigma = %.1f',StimLevels(cc)),'FontSize',fontSizeTitle);
    
    axis square
    cb = colorbar;
    caxis(all_fr_range);
    ylabel(cb,'spks/s');
    
    if(cc < 4)
        cb.Visible = 'off';
    end
    
    xlim(xrange)
    ylim(yrange)
    set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
end



set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*4 PlotSize*1],'PaperPosition',[0 0 PlotSize*4 PlotSize*1]);
if(saveFigs)
    saveas(gcf,sprintf('%s/HH_firingRates.pdf', figDir));
end

for tt = 1:length(tts)
    figure(10+tt);
    clf
    for hh = 1:length(hhs)
        for cc = 1:4
            subplot(length(hhs),4,cc+(hh-1)*4);
            dotPlot(Gs,Gs,fits_gain.glm.r_bar(:,:,cc,tts(tt),hhs(hh))',[],glm_fr_range);

            if(cc == 1)
                xlabel('G_K [pS/{\mu}m^2]','FontSize',fontSizeLabel);
                ylabel('G_{Na} [pS/{\mu}m^2]','FontSize',fontSizeLabel);
            end
            if(hh == 1 && cc == 1)
                title(sprintf('stimulus: \\sigma = %.1f\n trained on %s\n hist filt: %dms',StimLevels(cc),TrainStrings{tts(tt)},spkHistLengths(hhs(hh))),'FontSize',fontSizeTitle);
            elseif(cc == 1)
                title(sprintf('stimulus: \\sigma = %.1f\n hist filt: %dms',StimLevels(cc),spkHistLengths(hhs(hh))),'FontSize',fontSizeTitle);
            else
                title(sprintf('stimulus: \\sigma = %.1f',StimLevels(cc)),'FontSize',fontSizeTitle);
            end

            axis square
            cb = colorbar;
            ylabel(cb,'spks/s');
            caxis(glm_fr_range);

            if(cc < 4)
                cb.Visible = 'off';
            end

            xlim(xrange)
            ylim(yrange)
            set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
        end
    end
    set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*4 PlotSize*length(hhs)],'PaperPosition',[0 0 PlotSize*4 PlotSize*length(hhs)]);
    if(saveFigs)
        saveas(gcf,sprintf('%s/GLM_firingRates_train%d.pdf', figDir,tts(tt)));
    end
end


%%




Gs = 500 + (1:15)*100;

ratios  = (Gs./Gs');

ratios(isnan(gainDistance.hh.JS(:,:,1))) = nan;

[ratios,order] = sort(ratios(:));



figure(20);
clf;
for tt = 1:length(tts)
    for cc = 1:4

    end





