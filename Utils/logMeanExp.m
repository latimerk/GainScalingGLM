% Computes the log mean of exponentiated values.
%
%  function [y] = logMeanExp(x,dir)
%
%   x   = log values to compute over
%   dir = dimension over which mean is computed (default 1 for matrices)
%
%   y   = computationally safer log(mean(exp(x,dir))
%
function [y] = logMeanExp(x,dir)
if(nargin < 2)
    dir = 1;
    if(size(x,1) == 1)
        dir = 2;
    end
end
y = logSumExp(x,dir) - log(size(x,dir));