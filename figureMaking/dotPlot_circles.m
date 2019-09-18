function [] = dotPlot_circles(xs,ys,colors,lw)

if(nargin < 4)
    lw = 2;
end



rr = min(min(diff(sort(ys))),min(diff(sort(xs))))*0.49; 
            
hold on
for ii = 1:length(xs)
    for jj = 1:length(ys)
        cc = colors(:,ii,jj)';
        if(sum(isnan(cc))==0)
            
            
            
            rectangle('Position',[-rr+xs(ii) -rr+ys(jj) 2*rr 2*rr],'LineStyle','-','EdgeColor',cc,'FaceColor','none','Curvature',[1 1],'LineWidth',lw);
        end
    end
end
hold off