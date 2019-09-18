% Computes the log sum of exponentiated values.
%
%  function [y] = logSumExp(x,dir)
%
%   x   = log values to compute over
%   dir = dimension over which sum is computed (default 1 for matrices)
%
%   y   = computationally safer log(sum(exp(x,dir))
%
function [y] = logSumExp(x,dir)

if(nargin < 2)
    dir = 1;
    if(size(x,1) == 1)
        dir = 2;
    end
end

c = max(x,[],dir);
        
y = log(nansum(exp(x-c),dir)) + c;