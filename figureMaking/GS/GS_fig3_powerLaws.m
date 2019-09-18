load Results/Mease/Mease_gainInfo.mat
load Results/Mease/MeaseGLMs.mat
load('/home/latimerk/gitCode/GainScaling/Results/Mease/gainDistance.mat')
load('/home/latimerk/gitCode/GainScaling/Results/Mease/glmFitMetrics.mat')
figDir = 'Figs/GainScaling/';

addpath ~/FigureComposer/
addpath ~/FigureComposer/matlab/

addpath figureMaking/

saveFigs = false;


fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
PlotSize = 3;
         
spkHistLengths(1) = 10;

H = 16;
D = 3;

gg_c = nan(15,15,5);
for gg = 1:5
    gd = load(sprintf('/home/latimerk/gitCode/GainScaling/Results/Mease/GLMs_p%d/gainDistance.mat',gg*10),'gainDistance');
    gg_c(:,:,gg) = gd.gainDistance.glm.EM(:,:,D,H,T);
end

%%

Gs = 600:100:2000;
ratios  = (Gs./Gs')';
ratios(isnan(gainDistance.hh.EM(:,:,1))) = nan;
[ratios,order] = sort(ratios(:));

T = 5;

figure(1);
clf;

NR = 1;
NC = 4;


pps = 2:5;
cs = repmat(linspace(0.8,0.0,length(pps))',[1 3]);


subplot(NR,NC,1)
xx = -1:0.01:4;

semilogy(xx,exp(xx),'r','LineWidth',1);
hold on


for pp = 1:length(pps)
    if(pps(pp) < 6)
        plot(xx,log1p(exp(xx)).^pps(pp),'linewidth',1,'Color',cs(pp,:));
    end
end
xlim([xx(1),xx(end)])
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off



subplot(NR,NC,2)
hold on

tc = 1;

bps_0 = squeeze(mean(pseudo_R2_test(tc,:,:,T,H,1),1));

% plot(ratios,bps_0(order),'b');
plot([0 2],[0 0],'r:');
for pp = 1:length(pps)
    bps = squeeze(mean(pseudo_R2_test(tc,:,:,T,H,pps(pp)+1),1))-bps_0;
    plot(ratios,bps(order),'Color',cs(pp,:));
end
ylim([-0.1 0.1]);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off



subplot(NR,NC,3)
hold on

tc = 4;

bps_0 = squeeze(mean(pseudo_R2_test(tc,:,:,T,H,1),1));

% plot(ratios,bps_0(order),'b');
plot([0 2],[0 0],'r:');
for pp = 1:length(pps)
    bps = squeeze(mean(pseudo_R2_test(tc,:,:,T,H,pps(pp)+1),1))-bps_0;
    plot(ratios,bps(order),'Color',cs(pp,:));
end
ylim([-0.1 0.1]);
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off


subplot(NR,NC,4)
hold on
gg_hh = gainDistance.hh.EM(:,:,D);
gg_0 = gainDistance.glm.EM(:,:,D,H,T);


plot(ratios,gg_hh(order),'b');
plot(ratios,gg_0(order),'r');

for pp = 1:length(pps)
    gg_cp = gg_c(:,:,pps(pp));
    plot(ratios,gg_cp(order),'Color',cs(pp,:));
end
set(gca,'tickdir','out','box','off','fontsize',fontSizeAxis);
hold off


set(gcf,'PaperSize',[NC NR]*PlotSize,'PaperPosition',[0 0 NC NR]*PlotSize);

if(saveFigs)
    matfig2fyp(gcf,sprintf('%s/fig3_powerLaw_raw.fyp', figDir));
end