alpha = 0;


tts = 5;
hhs = [1 4 8 11 16];

if(alpha <= 0)
    saveDir = sprintf('Results/Mease/');
else
    saveDir = sprintf('Results/Mease/GLMs_p%d/',floor(alpha*10));
end

load(sprintf('%s/Mease_gainInfo.mat',saveDir));
load(sprintf('%s/MeaseGLMs.mat',saveDir));

load('Results/Mease/MeaseSims.mat', 'stim_mu','xs','frFunctions', 'G_mat_info_paper')

spkHistLengths(1) = 10;

fr_0 = squeeze(frFunctions(xs==0,:,:));

addpath figureMaking/

saveFigs = false;

%%
distanceFile = sprintf('%s/gainDistance.mat',saveDir);

fprintf('loading distance metrics...\n');
load(distanceFile);


%%

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 3;

Gs = 500 + (1:15)*100;

xticks = 600:400:2000;
yticks = 600:400:2000;

xrange = [Gs(1)-50 Gs(end)+50];
yrange = [Gs(1)-50 Gs(end)+50];


hh_JS_range = [0 ceil(nanmax(gainDistance.hh.JS(:))*100)/100];
hh_EM_range = [0 ceil(nanmax(gainDistance.hh.EM(:))*100)/100];


glm_JS_range = [0 ceil(nanmax(reshape(gainDistance.glm.JS(:,:,:,hhs,tts),[],1))*100)/100];
glm_EM_range = [0 ceil(nanmax(reshape(gainDistance.glm.EM(:,:,:,hhs,tts),[],1))*100)/100];

all_EM_range = [0 max(hh_EM_range(2),glm_EM_range(2))];
all_JS_range = [0 max(hh_JS_range(2),glm_JS_range(2))];

TrainStrings = cell(length(StimLevels)+1,1);
for ii = 1:length(StimLevels)
    TrainStrings{ii} = sprintf('\\sigma = %.1f',StimLevels(ii));
end
TrainStrings{end} = 'all \sigma';

if(alpha <= 0)
    figDir = sprintf('Figs/GainScaling/Distances/');
else
    figDir = sprintf('Figs/GainScaling/Distances/GLMs_p%d/',floor(alpha*10));
end
if(~isfolder(figDir))
    mkdir(figDir);
end
%%

figure(1);
clf
for cc = 1:3
    subplot(2,3,cc);
    dotPlot(Gs,Gs,gainDistance.hh.JS(:,:,cc)',[],all_JS_range);

    if(cc == 1)
        xlabel('G_K [pS/{\mu}m^2]','FontSize',fontSizeLabel);
        ylabel('G_{Na} [pS/{\mu}m^2]','FontSize',fontSizeLabel);
    end
    title(sprintf('JS: \\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);

    axis square
    cb = colorbar;
    ylabel(cb,'D_{\sigma} (bits)');
    caxis(all_JS_range);

    if(cc < 3)
        cb.Visible = 'off';
    end

    xlim(xrange)
    ylim(yrange)
    set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
end

for cc = 1:3
    subplot(2,3,3+cc);
    dotPlot(Gs,Gs,gainDistance.hh.EM(:,:,cc)',[],all_EM_range);

    if(cc == 1)
        xlabel('G_K [pS/{\mu}m^2]','FontSize',fontSizeLabel);
        ylabel('G_{Na} [pS/{\mu}m^2]','FontSize',fontSizeLabel);
    end
    title(sprintf('Wasserstein: \\sigma = %.1f',StimLevels(cc+1)),'FontSize',fontSizeTitle);

    axis square
    cb = colorbar;
    ylabel(cb,'Wass.');
    caxis(all_EM_range);

    if(cc < 3)
        cb.Visible = 'off';
    end

    xlim(xrange)
    ylim(yrange)
    set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
end

set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*2],'PaperPosition',[0 0 PlotSize*3 PlotSize*2]);
if(alpha < 1 && saveFigs)
    saveas(gcf,sprintf('%s/HH_gainScaling.pdf', figDir));
end

%%

for tt = 1:length(tts)
    
    figure(10+tt);
    clf
    for hh = 1:length(hhs)
        for cc = 1:3
            subplot(length(hhs),3,cc+3*(hh-1));
            dotPlot(Gs,Gs,gainDistance.glm.JS(:,:,cc,hhs(hh),tts(tt))',[],all_JS_range);

            if(cc == 1)
                xlabel('G_K [pS/{\mu}m^2]','FontSize',fontSizeLabel);
                ylabel('G_{Na} [pS/{\mu}m^2]','FontSize',fontSizeLabel);
            end
            if(hh == 1)
                if(cc == 1)
                    title(sprintf('\\alpha = %.1f, JS: \\sigma = %.1f,\n train on %s\n hist filt: %dms',alpha,StimLevels(cc+1),TrainStrings{tts(tt)},spkHistLengths(hhs(hh))),'FontSize',fontSizeTitle);
                else
                    title(sprintf('JS: \\sigma = %.1f,\n train on %s',StimLevels(cc+1),TrainStrings{tts(tt)}),'FontSize',fontSizeTitle);
                end
            elseif(cc == 1)
                title(sprintf('hist filt: %dms',spkHistLengths(hhs(hh))),'FontSize',fontSizeTitle);
            end

            axis square
            cb = colorbar;
            ylabel(cb,'D_{\sigma} (bits)');
            caxis(all_JS_range);

            if(cc < 3 || hh > 1)
                cb.Visible = 'off';
            end

            xlim(xrange)
            ylim(yrange)
            set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
        end
    end
    drawnow;
    set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*length(hhs)],'PaperPosition',[0 0 PlotSize*3 PlotSize*length(hhs)]);
    if(saveFigs)
        saveas(gcf,sprintf('%s/GLM_gainScaling_JS_train%d.pdf', figDir,tts(tt)));
    end
end



for tt = 1:length(tts)
    figure(20+tt);
    clf
    for hh = 1:length(hhs)
        for cc = 1:3
            subplot(length(hhs),3,cc+3*(hh-1));
            dotPlot(Gs,Gs,gainDistance.glm.EM(:,:,cc,hhs(hh),tts(tt))',[],all_EM_range);

            if(cc == 1)
                xlabel('G_K [pS/{\mu}m^2]','FontSize',fontSizeLabel);
                ylabel('G_{Na} [pS/{\mu}m^2]','FontSize',fontSizeLabel);
            end
            if(hh == 1)
                if(cc == 1)
                    title(sprintf('\\alpha = %.1f, Wasserstein: \\sigma = %.1f,\n train on %s\n hist filt: %dms',alpha, StimLevels(cc+1),TrainStrings{tts(tt)},spkHistLengths(hhs(hh))),'FontSize',fontSizeTitle);
                else
                    title(sprintf('Wasserstein: \\sigma = %.1f,\n train on %s',StimLevels(cc+1),TrainStrings{tts(tt)}),'FontSize',fontSizeTitle);
                end
            elseif(cc == 1)
                title(sprintf('hist filt: %dms',spkHistLengths(hhs(hh))),'FontSize',fontSizeTitle);
            end

            axis square
            cb = colorbar;
            ylabel(cb,'Wass.');
            caxis(all_EM_range);

            if(cc < 3 || hh > 1)
                cb.Visible = 'off';
            end

            xlim(xrange)
            ylim(yrange)
            set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
        end
    end
    
    drawnow;
    set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*3 PlotSize*length(hhs)],'PaperPosition',[0 0 PlotSize*3 PlotSize*length(hhs)]);
    if(saveFigs)
        saveas(gcf,sprintf('%s/GLM_gainScaling_EMD_train%d.pdf', figDir,tts(tt)));
    end
end




