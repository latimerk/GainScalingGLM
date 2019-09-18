alphas = [0 1.0 2.0 3.0 4.0 5.0 ];

alphas_range = [min(alphas(alphas>0)) max(alphas)];

load('Results/Mease/MeaseSims.mat', 'stim_mu','xs','frFunctions', 'G_mat_info_paper')
load('Results/Mease/MeaseGLMs.mat');
spkHistLengths(1) = 10;

tt = 5;

hhs = [4 8 12 16];

NA = length(alphas);
HH = length(hhs);

LLs  = nan(15,15,NA,HH);

for aa = 1:NA
    alpha = alphas(aa);
    
    if(alpha <= 0)
        saveDir = sprintf('Results/Mease/');
    else
        saveDir = sprintf('Results/Mease/GLMs_p%d/',floor(alpha*10));
    end

    %load(sprintf('%s/Mease_gainInfo.mat',saveDir));
    load(sprintf('%s/MeaseGLMs.mat',saveDir));
    
    glm_ll_tot = squeeze(sum(glm_ll(:,:,:,tt,hhs),1));
    LLs(:,:,aa,:) = glm_ll_tot;
end


for nn = 1:15
    for kk = 1:15
        if(frFunctions(xs == 0,nn,kk) > 1)
            LLs(nn,kk,:,:) = nan;
        end
    end
end

%%

figure(1);
clf;


fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 3;

Gs = 600:100:2000;

xticks = 600:400:2000;
yticks = 600:400:2000;

xrange = [Gs(1)-50 Gs(end)+50];
yrange = [Gs(1)-50 Gs(end)+50];

for hh = 1:HH
    subplot(1,HH,hh);
    hold on
    [~,mm] = max(LLs(:,:,:,hh),[],3);
    bb = nan(15,15);
    for ii = 1:15
        for jj = 1:15
            bb(ii,jj) = alphas(mm(ii,jj));
        end
    end
    bb(isnan(LLs(:,:,1,hh))) = nan;
    
    dotPlot(Gs,Gs,bb',[],alphas_range,[0 0 0]);
    
    title(sprintf('spk hist length: %dms',spkHistLengths(hhs(hh))),'FontSize',fontSizeTitle);

    axis square
    cb = colorbar;
    ylabel(cb,'best GLM nonlin.');
    caxis(alphas_range);

    if(hh < HH)
        cb.Visible = 'off';
    end

    xlim(xrange)
    ylim(yrange)
    set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
end