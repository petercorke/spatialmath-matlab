%PLOT_CIRCLE Draw a circle
%
% plot_circleC, R, OPTIONS) draws a circle on the current plot with 
% centre C=[X,Y] and radius R.  If C=[X,Y,Z] the circle is drawn in the
% XY-plane at height Z.
%
% If C (2xN) then N circles are drawn.  If R (1x1) then all
% circles have the same radius or else R (1xN) to specify the radius of
% each circle.
%
% H = plot_circle(...) as above but return handles. For multiple
% circles H is a vector of handles, one per circle.
%
% Options::
% 'edgecolor'   the color of the circle's edge, Matlab color spec
% 'fillcolor'   the color of the circle's interior, Matlab color spec
% 'alpha'       transparency of the filled circle: 0=transparent, 1=solid
% 'alter',H     alter existing circles with handle H
%
% - For an unfilled circle:
%   - any standard MATLAB LineStyle such as 'r' or 'b---'.
%   - any MATLAB LineProperty options can be given such as 'LineWidth', 2.
% - For a filled circle any MATLAB PatchProperty options can be given.
%
% Example::
%
%          H = plot_circle([3 4]', 2, 'r');  % draw red circle
%          plot_circle([3 4]', 3, 'alter', H); % change the circle radius
%          plot_circle([3 4]', 3, 'alter', H, 'LineColor', 'k'); % change the color
%
% Notes::
% - The 'alter' option can be used to create a smooth animation.
% - The circle(s) is added to the current plot irrespective of hold status.
%
% See also PLOT_ELLIPSE, PLOT_BOX, PLOT_POLY.

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
function handles = plot_circle(centre, rad, varargin)

    opt.fillcolor = [];
    opt.alpha = 1;
    opt.edgecolor = 'k';
    opt.alter = [];

    [opt,arglist] = tb_optparse(opt, varargin);
    
    if ~isempty(opt.alter) & ~ishandle(opt.alter)
        error('SMTB:plot_circle:badarg', 'argument to alter must be a valid graphic object handle');
    end

    holdon = ishold;
    hold on

	n = 50;
	th = [0:n]'/n*2*pi;
    
    if length(rad) == 1
        rad = rad*ones(numcols(centre),1);
    end
    if length(centre) == 2 || length(centre) == 3
        centre = centre(:);
    end

    for i=1:numcols(centre)
        x = rad(i)*cos(th) + centre(1,i);
        y = rad(i)*sin(th) + centre(2,i);
        if numrows(centre) > 2
            % plot 3D data
            z = ones(size(x))*centre(3,i);
            if isempty(opt.alter)
                h(i) = plot3(x, y, z, varargin{:});
            else
                set(opt.alter(i), 'xdata', x, 'ydata', y, 'zdata', z, arglist{:});
            end
        else
            % plot 2D data
            if isempty(opt.fillcolor)
                if isempty(opt.alter)
                    h(i) = plot(x, y, arglist{:});
                else
                    set(opt.alter(i), 'xdata', x, 'ydata', y, arglist{:});
                end
            else
                if isempty(opt.alter)
                    h(i) = patch(x, y, 0*y, 'FaceColor', opt.fillcolor, ...
                        'FaceAlpha', opt.alpha, 'EdgeColor', opt.edgecolor, arglist{:});
                else
                    set(opt.alter(i), 'xdata', x, 'ydata', y, arglist{:});
                end
                
            end
        end
    end

    if holdon == 0
        hold off
    end
    
    if nargout > 0
        handles = h;
    end
