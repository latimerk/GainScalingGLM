function [r2,r,mse] = getR2(y_data,y_model)

if(numel(y_data) ~= numel(y_model))
    error('Array sizes do not match');
end

vv = ~isnan(y_data(:)) & ~isnan(y_model(:));
y_data = y_data(vv);
y_model = y_model(vv);

if(sum(vv) >= 5)

    [a,~] = corrcoef(y_data,y_model);
    r = a(1,2);

    r2 = 1 - sum((y_data-y_model).^2)/sum((y_data-mean(y_data)).^2);

    mse = mean(abs((y_data-y_model)));
    
else
    r2 = nan;
    r = nan;
    mse = nan;
end