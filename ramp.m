%RAMP   create a ramp vector
%
% RAMP(N) output a vector of length N that ramps linearly from 0 to 1
%
%
% RAMP(V) as above but vector is same length as V
%
% RAMP(V, D) as above but elements increment by D.
%
% See also LINSPACE.
function r = ramp(v, d)
    if isscalar(v),
        l = v;
    else
        l = length(v);
    end
    if nargin == 1,
        r = [0:(l-1)]'/(l-1);
    else
        r = [0:(l-1)]'*d;
    end
