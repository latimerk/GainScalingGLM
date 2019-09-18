load('Results/Mease/HHexample.mat');
figDir = 'Figs/GainScaling/';

saveFigs = false;
tts_idx = (100*downsampleRate)+(1:(500*downsampleRate));
tts = (1:length(tts_idx))./downsampleRate;

xlims = [0 tts(end)];

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 3;


figure(1);
clf

s2_color = [160 0 15]./255;

if(~exist('fits_gain','var'))
    load('/home/latimerk/gitCode/GainScaling/Results/Mease/Mease_gainInfo.mat', 'fits_gain');
end

Gs = 600:100:2000;

[~,G_Na_idx] = min(abs(Gs-G_Na));
[~,G_K_idx]  = min(abs(Gs-G_Na));

T_sta = 100;

sta_1 = fits_gain.hh.sta_n(1:T_sta,G_Na_idx,G_K_idx,1);
sta_2 = fits_gain.hh.sta_n(1:T_sta,G_Na_idx,G_K_idx,4);

fy1 = conv(y,sta_1);
fy2 = conv(y,sta_2);

%% sigma 1 voltage
subplot(4,4,1)
plot(tts,squeeze(V(1,1,tts_idx)),'k');

ylabel('mV','FontSize',fontSizeLabel);
xlabel('time (ms)','FontSize',fontSizeLabel);
ylim([-90 40]);
xlim(xlims);

set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
%% sigma 1 input

subplot(4,4,5);
plot(tts,y(tts_idx)*sigs(1) + stim_dc,'k');

title(sprintf('\\sigma = %.1f',StimLevels(1)),'FontSize',fontSizeTitle);
ylabel('\mu A/cm^2','FontSize',fontSizeLabel);
xlabel('time (ms)','FontSize',fontSizeLabel);

ylim([-6 6]);
xlim(xlims);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% sigma 1 voltage
subplot(4,4,1)
plot(tts,squeeze(V(1,1,tts_idx)),'k');

ylabel('mV','FontSize',fontSizeLabel);
xlabel('time (ms)','FontSize',fontSizeLabel);
title(sprintf('\\sigma = %.1f',StimLevels(1)),'FontSize',fontSizeTitle);
ylim([-90 40]);
xlim(xlims);

set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% sigma 2 voltage
subplot(4,4,2)
plot(tts,squeeze(V(1,4,tts_idx)),'Color',s2_color);

ylim([-90 40]);
title(sprintf('\\sigma = %.1f',StimLevels(4)),'FontSize',fontSizeTitle);
xlim(xlims);

set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% sigma 1 input

subplot(4,4,5);
plot(tts,y(tts_idx)*sigs(1) + stim_dc,'k');

ylabel('\mu A/cm^2','FontSize',fontSizeLabel);
xlabel('time (ms)','FontSize',fontSizeLabel);

ylim([-6 6]);
xlim(xlims);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);
%% sigma 2 input

subplot(4,4,6);
plot(tts,y(tts_idx)*sigs(4) + stim_dc,'Color',s2_color);


ylim([-6 6]);
xlim(xlims);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% stas

subplot(4,4,9:10);
hold on

plot(sta_1,'k');
plot(sta_2,'Color',s2_color);

hold off
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);


%% stas

subplot(4,4,9:10);
hold on

plot(sta_1,'k');
plot(sta_2,'Color',s2_color);

hold off
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

% %% filtered stim 1
% 
% subplot(3,4,7);
% plot(tts,fy1(tts_idx)*sigs(1),'k');
% 
% xlabel('time (ms)','FontSize',fontSizeLabel);
% 
% ylim([-10 10]);
% xlim(xlims);
% set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

% %% filtered stim 2
% 
% subplot(3,4,8);
% plot(tts,fy2(tts_idx)*sigs(4),'Color',s2_color);
% 
% xlabel('time (ms)','FontSize',fontSizeLabel);
% 
% ylim([-10 10]);
% xlim(xlims);
% set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% p_spk

subplot(4,4,7);
hold on
dx1 = (fits_gain.bins(2)-fits_gain.bins(1))*sigs(1);
dx2 = (fits_gain.bins(2)-fits_gain.bins(1))*sigs(4);

plot(fits_gain.bins*sigs(1),fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,1)./dx1,'Color',[0 0 0]);
plot(fits_gain.bins*sigs(4),fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,4)./dx2,'Color',s2_color);

prior1 = fits_gain.hh.s_all(:,G_Na_idx,G_K_idx,1);
prior2 = fits_gain.hh.s_all(:,G_Na_idx,G_K_idx,4);
prior1 = prior1(:)./(sum(prior1)*dx1);
prior2 = prior2(:)./(sum(prior2)*dx2);

plot(fits_gain.bins*sigs(1),prior1,'Color','b','LineWidth',0.5);
plot(fits_gain.bins*sigs(4),prior2,'Color','g','LineWidth',0.5);

xlabel('s','FontSize',fontSizeLabel);
ylabel('p(s|spk)','fontsize',fontSizeLabel);

xlim([-8 8]);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% p_spk scaled

subplot(4,4,8);
hold on
dx1 = (fits_gain.bins(2)-fits_gain.bins(1));
dx2 = (fits_gain.bins(2)-fits_gain.bins(1));

plot(fits_gain.bins,fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,1)./dx1,'Color',[0 0 0]);
plot(fits_gain.bins,fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,4)./dx2,'Color',s2_color);

prior1 = fits_gain.hh.s_all(:,G_Na_idx,G_K_idx,1);
prior2 = fits_gain.hh.s_all(:,G_Na_idx,G_K_idx,4);
prior1 = prior1(:)./(sum(prior1)*dx1);
prior2 = prior2(:)./(sum(prior2)*dx2);

plot(fits_gain.bins,prior1,'Color','b','LineWidth',0.5);
plot(fits_gain.bins,prior2,'Color','g','LineWidth',0.5);

xlabel('\hat{s}','FontSize',fontSizeLabel);
ylabel('p(\hat{s}|spk)','fontsize',fontSizeLabel);

js = getJS(fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,1),fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,4))./dx1;

title(sprintf('JS = %.2f',js))

xlim([-4 4]);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% R

subplot(4,4,11);
dx1 = (fits_gain.bins(2)-fits_gain.bins(1))*sigs(1);
dx2 = (fits_gain.bins(2)-fits_gain.bins(1))*sigs(4);


like1 = fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,1)./dx1;
like2 = fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,4)./dx2;
prior1 = fits_gain.hh.s_all(:,G_Na_idx,G_K_idx,1);
prior2 = fits_gain.hh.s_all(:,G_Na_idx,G_K_idx,4);
prior1 = prior1(:)./(sum(prior1)*dx1);
prior2 = prior2(:)./(sum(prior2)*dx2);

semilogy(fits_gain.bins,like1./prior1,'Color',[0 0 0]);
hold on
semilogy(fits_gain.bins*sigs(4),like2./prior2,'Color',s2_color);


xlabel('s','FontSize',fontSizeLabel);
ylabel('R_\sigma(s)','fontsize',fontSizeLabel);

xlim([-8 8]);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% R scaled

subplot(4,4,12);
dx1 = (fits_gain.bins(2)-fits_gain.bins(1));
dx2 = (fits_gain.bins(2)-fits_gain.bins(1));


like1 = fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,1)./dx1;
like2 = fits_gain.hh.p_spk(:,G_Na_idx,G_K_idx,4)./dx2;
prior1 = fits_gain.hh.s_all(:,G_Na_idx,G_K_idx,1);
prior2 = fits_gain.hh.s_all(:,G_Na_idx,G_K_idx,4);
prior1 = prior1(:)./(sum(prior1)*dx1);
prior2 = prior2(:)./(sum(prior2)*dx2);

semilogy(fits_gain.bins,like1./prior1,'Color',[0 0 0]);
hold on
semilogy(fits_gain.bins,like2./prior2,'Color',s2_color);


xlabel('\hat{s}','FontSize',fontSizeLabel);
ylabel('R_\sigma(\hat{s})','fontsize',fontSizeLabel);

xlim([-4 4]);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

G_Na_idx2 = 1;
G_K_idx2  = 15;
sta_1 = fits_gain.hh.sta_n(1:T_sta,G_Na_idx2,G_K_idx2,1);
sta_2 = fits_gain.hh.sta_n(1:T_sta,G_Na_idx2,G_K_idx2,4);

%% p_spk scaled 2

subplot(4,4,15);
hold on
dx1 = (fits_gain.bins(2)-fits_gain.bins(1));
dx2 = (fits_gain.bins(2)-fits_gain.bins(1));

plot(fits_gain.bins,fits_gain.hh.p_spk(:,G_Na_idx2,G_K_idx2,1)./dx1,'Color',[0 0 0]);
plot(fits_gain.bins,fits_gain.hh.p_spk(:,G_Na_idx2,G_K_idx2,4)./dx2,'Color',s2_color);

prior1 = fits_gain.hh.s_all(:,G_Na_idx2,G_K_idx2,1);
prior2 = fits_gain.hh.s_all(:,G_Na_idx2,G_K_idx2,4);
prior1 = prior1(:)./(sum(prior1)*dx1);
prior2 = prior2(:)./(sum(prior2)*dx2);

plot(fits_gain.bins,prior1,'Color','b','LineWidth',0.5);
plot(fits_gain.bins,prior2,'Color','g','LineWidth',0.5);

xlabel('\hat{s}','FontSize',fontSizeLabel);
ylabel('p(\hat{s}|spk)','fontsize',fontSizeLabel);
js = getJS(fits_gain.hh.p_spk(:,G_Na_idx2,G_K_idx2,1),fits_gain.hh.p_spk(:,G_Na_idx2,G_K_idx2,4))./dx1;

title(sprintf('JS = %.2f',js))

xlim([-4 4]);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);

%% R scaled 2

subplot(4,4,16);
dx1 = (fits_gain.bins(2)-fits_gain.bins(1));
dx2 = (fits_gain.bins(2)-fits_gain.bins(1));


like1 = fits_gain.hh.p_spk(:,G_Na_idx2,G_K_idx2,1)./dx1;
like2 = fits_gain.hh.p_spk(:,G_Na_idx2,G_K_idx2,4)./dx2;
prior1 = fits_gain.hh.s_all(:,G_Na_idx2,G_K_idx2,1);
prior2 = fits_gain.hh.s_all(:,G_Na_idx2,G_K_idx2,4);
prior1 = prior1(:)./(sum(prior1)*dx1);
prior2 = prior2(:)./(sum(prior2)*dx2);

semilogy(fits_gain.bins,like1./prior1,'Color',[0 0 0]);
hold on
semilogy(fits_gain.bins,like2./prior2,'Color',s2_color);


xlabel('\hat{s}','FontSize',fontSizeLabel);
ylabel('R_\sigma(\hat{s})','fontsize',fontSizeLabel);

xlim([-4 4]);
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);


%% stas

subplot(4,4,13:14);
hold on

plot(sta_1,'k');
plot(sta_2,'Color',s2_color);

hold off
set(gca,'TickDir','out','box','off','XTick',xticks,'YTick',yticks,'FontSize',fontSizeAxis);


%%

set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*4 PlotSize*2],'PaperPosition',[0 0 PlotSize*4 PlotSize*2]);
if(saveFigs)
    matfig2fyp(gcf,sprintf('%s/gainScaling_diagram.fyp',figDir));
    %saveas(gcf,sprintf('%s/exampleHH.pdf', figDir));
end