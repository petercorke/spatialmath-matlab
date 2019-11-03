%ISVEC Test if vector
%
% ISVEC(V) is true (1) if the argument V is a 3-vector, either a 
% row- or column-vector.  Otherwise false (0).
%
% ISVEC(V, L) is true (1) if the argument V is a vector of length L,
% either a row- or column-vector.  Otherwise false (0).
%
% Notes::
% - Differs from MATLAB builtin function ISVECTOR which returns true
%   for the case of a scalar, ISVEC does not.
% - Gives same result for row- or column-vector, ie. 3x1 or 1x3 gives true.
% - Works for a symbolic math symfun.
%
% See also ISHOMOG, ISROT.

% Copyright (C) 1993-2019 Peter I. Corke
%
% This file is part of The Spatial Math Toolbox for MATLAB (SMTB).
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
% of the Software, and to permit persons to whom the Software is furnished to do
% so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
% IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% https://github.com/petercorke/spatial-math

function h = isvec(v, l)
    if nargin == 1
            l = 3;
    end
    if isa(v, 'symfun')
        h = logical( length(formula(v)) == l);
    else
        d = size(v);
        h = logical( length(d) == 2 && min(d) == 1 && numel(v) == l );
    end
end

