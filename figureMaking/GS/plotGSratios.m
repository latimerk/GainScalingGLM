alpha = 3;
if(alpha <= 0)
    saveDir = sprintf('Results/Mease/');
else
    saveDir = sprintf('Results/Mease/GLMs_p%d/',floor(alpha*10));
end

load(sprintf('%s/gainDistance.mat',saveDir));
load(sprintf('%s/MeaseGLMs.mat',saveDir));

spkHistLengths(1) = 10;

saveFigs = false;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 3;

Gs = 500 + (1:15)*100;

ratios  = (Gs./Gs')';

ratios(isnan(gainDistance.hh.JS(:,:,1))) = nan;


xx = [0.2 2.0];
xt = 0.2:0.6:2.0;

% xx = [-0.3 0.6];
% xt = -0.3:0.3:6;


[ratios,order] = sort(ratios(:));


tts = [1 4 5];
hhs = [1  3 4 6 11 16];


hh_JS_range = [0 ceil(nanmax(gainDistance.hh.JS(:))*100)/100];
hh_EM_range = [0 ceil(nanmax(gainDistance.hh.EM(:))*100)/100];


glm_JS_range = [0 ceil(nanmax(reshape(gainDistance.glm.JS(:,:,:,hhs,tts),[],1))*100)/100];
glm_EM_range = [0 ceil(nanmax(reshape(gainDistance.glm.EM(:,:,:,hhs,tts),[],1))*100)/100];

all_EM_range = [0 max(hh_EM_range(2),glm_EM_range(2))];
all_JS_range = [0 0.07];%max(hh_JS_range(2),glm_JS_range(2))];

lineWidth_hh = 2;
lineStyle_hh = '-';
lineColor_hh = [0 0 1;
                0 0 1;
                0 0 1];
            
            
lineWidth_glm = 1;
lineStyle_glm = '-';
lineColor_glm = ones(length(hhs),3).*linspace(0.8,0,length(hhs))';

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
for cc = 3
    figure(cc);
    clf;
    for tt = 1:length(tts)
        subplot(2,length(tts),tt);
        hold on

        gg = gainDistance.hh.JS(:,:,cc);
        gg = gg(order);

        plot(ratios,gg,lineStyle_hh,'LineWidth',lineWidth_hh,'Color',lineColor_hh(cc,:));


        for hh = 1:length(hhs)
            gg = gainDistance.glm.JS(:,:,cc,hhs(hh),tts(tt));
            gg = gg(order);

            plot(ratios,gg,lineStyle_glm,'LineWidth',lineWidth_glm,'Color',lineColor_glm(hh,:));
        end

        if(tt == 1)
            title(sprintf('\\alpha=%.1f, \\sigma = %.1f\n train on %s',alpha,StimLevels(cc+1),TrainStrings{tts(tt)}),'FontSize',fontSizeTitle);
            xlabel('G_{Na}/G_{K} ratio','FontSize',fontSizeLabel);
            ylabel('JS divergence','FontSize',fontSizeLabel);
        else
            title(sprintf('train on %s',TrainStrings{tts(tt)}),'FontSize',fontSizeTitle);
        end

        if(tt == 1)
            legendObjs = cell(length(hhs)+1,1);
            legendObjs{1} = 'Hodgkin-Huxley';
            for hh = 1:length(hhs)
                legendObjs{hh+1} = sprintf('GLM: spk hist %d ms',spkHistLengths(hhs(hh)));
            end
            legend(legendObjs,'Location','NorthWest','fontSize',fontSizeLegend);
        end
        
        hold off
        xlim(xx);
        ylim(all_JS_range);
        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis,'XTick',xt);     
    end

    for tt = 1:length(tts)
        subplot(2,length(tts),tt+length(tts));

        gg = gainDistance.hh.EM(:,:,cc);
        gg = gg(order);

        plot(ratios,gg,lineStyle_hh,'LineWidth',lineWidth_hh,'Color',lineColor_hh(cc,:));
        hold on


        for hh = 1:length(hhs)
            gg = gainDistance.glm.EM(:,:,cc,hhs(hh),tts(tt));
            gg = gg(order);

            plot(ratios,gg,lineStyle_glm,'LineWidth',lineWidth_glm,'Color',lineColor_glm(hh,:));
        end

        if(tt == 1)
            title(sprintf('\\sigma = %.1f\n train on %s',StimLevels(cc+1),TrainStrings{tts(tt)}),'FontSize',fontSizeTitle);
            xlabel('G_{Na}/G_{K} ratio','FontSize',fontSizeLabel);
            ylabel('Wasserstein distance','FontSize',fontSizeLabel);
        else
            title(sprintf('train on %s',TrainStrings{tts(tt)}),'FontSize',fontSizeTitle);
        end

        hold off
        xlim(xx);
        ylim(all_EM_range);
        set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis,'XTick',xt);     
    end
    
    set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*2 PlotSize*2],'PaperPosition',[0 0 PlotSize*2 PlotSize*2]);
    if(saveFigs)
        saveas(gcf,sprintf('%s/conductanceRatios_%d.pdf', figDir,cc));
    end
end



%%