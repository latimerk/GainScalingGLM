function [cm,range] = dotPlot2(xs,ys,zs,cm,range)


if(nargin < 4 || isempty(cm))
    cm = colormap;
end
if(nargin < 5)
    range = [nanmin(zs(:)),nanmax(zs(:))];
end

ss = linspace(range(1),range(2),size(cm,1));
ss(1) = -inf;
ss(end) = inf;

hold on
rr = min(min(diff(sort(ys))),min(diff(sort(xs))))*0.49; 
rec1 = rectangle('Position',[0 0 0 0],'LineStyle','none','FaceColor',[0 0 0],'Curvature',[1 1]);        
for ii = 1:length(xs)
    for jj = 1:length(ys)
        v = zs(ii,jj);
        if(~isnan(v))
            hh = histc(v,ss);
            cc = cm(hh>0,:);
            
            
            rec2 = copyobj(rec1,rec1.Parent);
            rec2.Position = [-rr+xs(ii) -rr+ys(jj) 2*rr 2*rr];
            rec2.FaceColor = cc;
        end
    end
end
hold off