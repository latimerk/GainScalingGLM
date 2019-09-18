alphas = [0 0.2 0.4 0.6];
NA = length(alphas);

xx1 = [zeros(50,1);ones(500,1);zeros(500,1);ones(500,1);zeros(400,1)];
T1 = 1;
R  = 2;
R2 = 4;
L = 32e3;
Fs = 1e3;
xx2 = sin((1:L)./(Fs)*2*pi./T1)';
xx3 = sin((1:L)./(Fs)*2*pi./(T1*R))';
xx4 = sin((1:L)./(Fs)*2*pi./(T1*R2))';

FF = Fs*(0:(L/2))/L;
FF = [FF FF(end-1:-1:2)]';
ww = FF*2*pi;

fontSizeTitle = 12;
fontSizeLabel = 12;
fontSizeAxis  = 10;
fontSizeLegend= 10;
PlotSize = 1;

NC = 3;

addpath Utils/

figure(1);
clf;

for aa = 1:NA
    subplot(NA,NC, 1 + (aa-1)*NC);
    hold on
    yy1 = fgl_deriv( alphas(aa), xx1, 1e-3);
    plot(yy1,'k');
    hold off
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);

    title(sprintf('alpha = %.1f',alphas(aa)));
    
    subplot(NA,NC, 2 + (aa-1)*NC);
    hold on
%     yy2 = fgl_deriv( alphas(aa), xx2, 1e-3);
    ff2 = fft(xx2);
    angle_new = (angle(ff2)+pi/2*alphas(aa));
    angle_new(L/2+2:end) = (angle(ff2(L/2+2:end)-pi/2*alphas(aa)));
    mag_new = abs(ff2).*ww.^alphas(aa);
    ff_new = (cos(angle_new) + sin(angle_new)*1i).*mag_new;
    yy2 = real(ifft(ff_new));
    yy2 = yy2(T1*Fs+(1:(T1*Fs)));
    
    plot(yy2,'k');
%     yy3 = fgl_deriv( alphas(aa), xx3, 1e-3);
%     yy3 = yy3(T1*R*Fs+(1:(T1*R*Fs)));
    ff3 = fft(xx3);
    angle_new = (angle(ff3)+pi/2*alphas(aa));
    angle_new(L/2+2:end) = (angle(ff3(L/2+2:end)-pi/2*alphas(aa)));
    mag_new = abs(ff3).*ww.^alphas(aa);
    ff_new = (cos(angle_new) + sin(angle_new)*1i).*mag_new;
    yy3 = real(ifft(ff_new));
    yy3 = yy3(T1*R*Fs+(1:(T1*R*Fs)));
    
    plot(yy3,'b');
%     yy4 = fgl_deriv( alphas(aa), xx4, 1e-3);
    ff4 = fft(xx4);
    angle_new = (angle(ff4)+pi/2*alphas(aa));
    angle_new(L/2+2:end) = (angle(ff4(L/2+2:end)-pi/2*alphas(aa)));
    mag_new = abs(ff4).*ww.^alphas(aa);
    ff_new = (cos(angle_new) + sin(angle_new)*1i).*mag_new;
    yy4 = real(ifft(ff_new));
    yy4 = yy4(T1*R2*Fs+(1:(T1*R2*Fs)));

    plot(yy4,'r');
    hold off
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);
    
    
    subplot(NA,NC, 3 + (aa-1)*NC);
    hold on
    plot(yy2,'k');
    
    yy3b = zeros(length(yy2),1);
    for ii = 1:length(yy2)
        yy3b(ii) = mean(yy3((ii-1)*R + (1:R)));
    end
    yy4b = zeros(length(yy2),1);
    for ii = 1:length(yy2)
        yy4b(ii) = mean(yy4((ii-1)*R2 + (1:R2)));
    end
    plot(yy3b,'b');
    plot(yy4b,'r');
    hold off
    set(gca,'TickDir','out','box','off','FontSize',fontSizeAxis);


end

set(gcf,'PaperUnits','inches','PaperSize',[PlotSize*NA PlotSize*NC],'PaperPosition',[0 0 PlotSize*NA PlotSize*NC]);
    