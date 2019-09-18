GS_fig1_example;
addpath figureMaking/
%%
PlotSize  = 3;

D = 3;
hh_ds  = gainDistance.hh.EM(:,:,D);
glm_ds = gainDistance.glm.EM(:,:,D,H,T);
c_range = [0 ceil(max(max(hh_ds(:)),max(glm_ds(:)))*20)/20];


c_range1 = [0 ceil(max(glm_ds(:))*20)/20];
c_range2 = [0 ceil(max(hh_ds(:))*20)/20];

lims = [600 2000]+[-1 1]*60;

exampleColors = [214 0 255;
                 230 0 0;
                 0   230 40]./255;
exampleColors(1,:) = 0;      
exampleColors(2,:) = 0;          
exampleColors(3,:) = 0;                 

dotColors = nan(3,15,15);
for ii = 1:NE
    dotColors(:,example_idx(ii,2),example_idx(ii,1)) = exampleColors(ii,:);
end



figure(3);
clf


subplot(2,1,1)
hold on
dotPlot(Gs,Gs,glm_ds',[],c_range);
dotPlot_circles(Gs,Gs,dotColors,3);
xlim(lims);
ylim(lims);
xlabel('G_{Na} (Ps/\mu m^2)','FontSize',fontSizeLabel);
ylabel('G_{K} (Ps/\mu m^2)','FontSize',fontSizeLabel);

cb1 = colorbar();
ylabel(cb1,'D_{\sigma=2}','FontSize',fontSizeLabel);
cb1.Visible = true;


set(gca,'XTick',600:400:2000,'YTick',600:400:2000);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off


subplot(2,1,2)
hold on
dotPlot(Gs,Gs,hh_ds',[],c_range);
dotPlot_circles(Gs,Gs,dotColors,3);
xlim(lims);
ylim(lims);
cb2 = colorbar();
ylabel(cb2,'D_{\sigma=2}','FontSize',fontSizeLabel);
cb2.Visible = false;

xlabel('G_{Na} (Ps/\mu m^2)','FontSize',fontSizeLabel);
ylabel('G_{K} (Ps/\mu m^2)','FontSize',fontSizeLabel);

set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);

set(gca,'XTick',600:400:2000,'YTick',600:400:2000);
hold off


%%
set(gcf,'PaperSize',[1 2]*PlotSize,'PaperPosition',[0 0 1 2]*PlotSize);

if(saveFigs)
    %matfig2fyp(gcf,sprintf('%s/fig1_GSsummary_raw.fyp', figDir));
    saveas(gcf,sprintf('%s/fig1_GSsummary_raw.eps', figDir),'epsc');
end