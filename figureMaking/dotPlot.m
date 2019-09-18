function [cm,range] = dotPlot(xs,ys,zs,cm,range,outOfRangeColor)


if(nargin < 4 || isempty(cm))
    cm = colormap;
end
if(nargin < 5)
    nn = zs(~isinf(zs(:)));
    range = [nanmin(nn),nanmax(nn)];
end
if(nargin < 6)
    outOfRangeColor = [cm(1,:); cm(end,:)];
end

ss = linspace(range(1),range(2),size(cm,1));
ss(1) = -inf;
ss(end) = inf;

rr = min(min(diff(sort(ys))),min(diff(sort(xs))))*0.49; 
            
hold on
for ii = 1:length(xs)
    for jj = 1:length(ys)
        v = zs(ii,jj);
        if(~isnan(v))
            if(v > range(2) || v < range(1))
                if(size(outOfRangeColor,1) > 1)
                    if(v > range(2))
                        cc = outOfRangeColor(2,:);
                    else
                        cc = outOfRangeColor(1,:);
                    end
                else
                    cc = outOfRangeColor(1,:);
                end
            else
                hh = histc(v,ss);
                cc = cm(hh>0,:);
            end
            
            rectangle('Position',[-rr+xs(ii) -rr+ys(jj) 2*rr 2*rr],'LineStyle','none','FaceColor',cc,'Curvature',[1 1]);
        end
    end
end
caxis(range);
hold off