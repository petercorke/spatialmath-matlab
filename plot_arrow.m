%PLOT_ARROW Draw an arrow in 2D or 3D
%
% PLOT_ARROW(P1, P2, OPTIONS) draws an arrow from P1 to P2 (2x1 or 3x1).  For 3D 
% case the arrow head is a cone.
%
% PLOT_ARROW(P, OPTIONS) as above where the columns of P (2x2 or 3x2) define the 
% start and end points, ie. P=[P1 P2].
%
% H = PLOT_ARROW(...) as above but returns the graphics handle of the arrow.
%
% Options::
% - All options are passed through to arrow3.  
% - MATLAB LineSpec such as 'r' or 'b--'
%
% Example::
%         plot_arrow([0 0 0]', [1 2 3]', 'r')  % a red arrow
%         plot_arrow([0 0 0], [1 2 3], 'r--3', 4) % dashed red arrow big head
%
% Notes::
% - Requires https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3
% - ARROW3 attempts to preserve the appearance of existing axes.  In
%   particular, ARROW3 will not change XYZLim, View, or CameraViewAngle.
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
function h = plot_arrow(p1, varargin)
    if isvec(p1, 2) || isvec(p1, 3)
        p1 = p1(:)';  % force row vector
    else
        error('SMTB:plot_arrow', 'P1 must have 2 or 3 elements')
    end

    if nargin > 1
        p2 = varargin{1};
        if isnumeric(p2)
            if isvec(p2, 2) || isvec(p2, 3)
                p2 = p2(:)';  % force row vector
                varargin = varargin(2:end);
             else
                error('SMTB:plot_arrow', 'P2 must have 2 or 3 elements')
            end
        end
    end
    assert(numcols(p1) == numcols(p2), 'SMTB:plot_arrow', 'P1 and P2 must be the same length');

    hh = arrow3(p1, p2, varargin{:});

    if nargout > 0
        h = hh;
    end
end
