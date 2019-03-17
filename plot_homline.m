%PLOT_HOMLINE Draw a line in homogeneous form
%
% PLOT_HOMLINE(L) draws a 2D line in the current plot defined in homogenous
% form ax + by + c = 0  where L (3x1) is L = [a b c].
% The current axis limits are used to determine the endpoints of
% the line.  If L (3xN) then N lines are drawn, one per column.
%
% PLOT_HOMLINE(L, LS) as above but the MATLAB line specification LS is given.  
%
% H = PLOT_HOMLINE(...) as above but returns a vector of graphics handles for the lines.
%
% Notes::
% - The line(s) is added to the current plot.
% - The line(s) can be drawn in 3D axes but will always lie in the
%   xy-plane.
%
% Example::
%          L = homline([1 2]', [3 1]'); % homog line from (1,2) to (3,1)
%          plot_homline(L, 'k--'); % plot dashed black line
%
% See also PLOT_BOX, PLOT_POLY, HOMLINE.

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


function handles = plot_homline(lines, varargin)

	% get plot limits from current graph
	xlim = get(gca, 'XLim');
	ylim = get(gca, 'YLim');

    ish = ishold;
    hold on;
    
    if min(size(lines)) == 1
        lines = lines(:);
    end
    
    assert(numrows(lines) == 3, 'SMTB:plot_homline:badarg', 'Input must be a 3-vector or 3xN matrix');

	h = [];
	% for all input lines (columns
	for l=lines
        if abs(l(2)) > abs(l(1))
            y = (-l(3) - l(1)*xlim) / l(2);
            hh = plot(xlim, y, varargin{:});
        else
            x = (-l(3) - l(2)*ylim) / l(1);
            hh = plot(x, ylim, varargin{:});
        end
		h = [h; hh];
	end

    if ~ish
        hold off
    end

	if nargout > 0,
		handles = h;
	end
