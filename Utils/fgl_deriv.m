function [Y,R] = fgl_deriv( a, y, h)% R )
% fgl_deriv
%
%   Computes the fractional derivative of order alpha (a) for the function
%   y sampled on a regular grid with spacing h, using the Grunwald-Letnikov
%   formulation.
%
%   Inputs:
%   a : fractional order
%   y : sampled function
%   h : period of the sampling lattice
%
%   Output:
%   Y : fractional derivative estimate given y
%
%   Note that this implementation is similar to that of Bayat 2007
%   (FileExchange: 13858-fractional-differentiator), and takes the exact
%   same inputs, but uses a vectorized formulation to accelerates the
%   computation in Matlab.
%
%   Copyright 2014 Jonathan Hadida
%   Contact: Jonathan dot hadida [a] dtc.ox.ac.uk

    n  = numel(y);
    J  = 0:(n-1);
    G1 = gamma( J+1 );
    G2 = gamma( a+1-J );
    s  = (-1) .^ J;

%     M  = tril( ones(n) );
%     if(nargin < 4)
%         R  = toeplitz( y,[y(1);zeros(n-1,1)]);
%     end
%     T  = meshgrid( (gamma(a+1)/(h^a)) * s ./ (G1.*G2) );
    T  =  ((gamma(a+1)/(h^a)) * s ./ (G1.*G2))' ;
%     Y  = reshape(sum( R .* T, 2 ), size(y));
%     Y  =  R * T;
    
    Y = conv(T,[zeros(n,1);y],'same');
    
    
end