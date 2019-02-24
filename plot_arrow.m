%PLOT_ARROW Draw an arrow in 2D or 3D
%
% PLOT_ARROW(P1, P2, OPTIONS) draws an arrow from P1 to P2 (2x1 or 3x1).
%
% PLOT_ARROW(P, OPTIONS) as above where the columns of P (2x2 or 3x2) define where P=[P1 P2].
%
% Options::
% - All options are passed through to arrow3.  
% - MATLAB colorspec such as 'r' or 'b--'
%
% See also ARROW3.

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
function plot_arrow(p1, varargin)
    
    if min(size(p1)) == 1
        % p1 is a vector
        p1 = p1(:);
        p2 = varargin{1};
        p2 = p2(:);
        assert(numrows(p1) == numrows(p2), 'SMTB:plot_arrow', 'P1 and P2 must be the same length');
        varargin = varargin{2:end};
    else
        % p1 is a 2-column matrix
        assert(numcols(p1) == 2, 'SMTB:plot_arrow', 'P1 must have 2 columns');
        p2 = p1(:,2);
        p1 = p1(:,1);
    end
    
    assert(any(numrows(p1) == [2 3]), 'SMTB:plot_arrow', '2D or 3D points only');

    arrow3(p1', p2', varargin{:});
