%PLOT_BOX Draw a box
%
% PLOT_BOX(B, OPTIONS) draws a box defined by B=[XL XR; YL YR] on the current
% plot with optional MATLAB linestyle options LS.
%
% PLOT_BOX(X1,Y1, X2,Y2, OPTIONS) draws a box with corners at (X1,Y1) and (X2,Y2),
% and optional MATLAB linestyle options LS.
%
% PLOT_BOX('centre', P, 'size', W, OPTIONS) draws a box with center at P=[X,Y] and
% with dimensions W=[WIDTH HEIGHT].
%
% PLOT_BOX('topleft', P, 'size', W, OPTIONS) draws a box with top-left at P=[X,Y] 
% and with dimensions W=[WIDTH HEIGHT].
%
% PLOT_BOX('matlab', BOX, LS) draws box(es) as defined using the MATLAB convention of
% specifying a region in terms of top-left coordinate, width and height.  One box is
% drawn for each row of BOX which is [xleft ytop width height].
%
%
% H = PLOT_ARROW(...) as above but returns the graphics handle of the arrow.
%
% Options::
% 'edgecolor'   the color of the circle's edge, MATLAB ColorSpec
% 'fillcolor'   the color of the circle's interior, MATLAB ColorSpec
% 'alpha'       transparency of the filled circle: 0=transparent, 1=solid
%
% - For an unfilled box:
%   - any standard MATLAB LineSpec such as 'r' or 'b---'.
%   - any MATLAB LineProperty options can be given such as 'LineWidth', 2.
% - For a filled box any MATLAB PatchProperty options can be given.
%
% Examples::
%          plot_box([0 1; 0 2], 'r')   % draw a hollow red box
%          plot_box([0 1; 0 2], 'fillcolor', 'b', 'alpha', 0.5)   % translucent filled blue box
%
% Notes::
% - The box is added to the current plot irrespective of hold status.
%
% See also PLOT_POLY, PLOT_CIRCLE, PLOT_ELLIPSE.

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

function h = plot_box(varargin)
    opt.centre = [];
    opt.topleft = [];
    opt.matlab = [];
    opt.size = [];
    opt.fillcolor = [];
    opt.alpha = 1;
    opt.edgecolor = 'k';

    [opt,args] = tb_optparse(opt, varargin);

    if ~isempty(opt.size)
        if size(opt.size) == 1
            w = opt.size;
            h = opt.size;
        else
            w = opt.size(1);
            h = opt.size(2);
        end

        if ~isempty(opt.centre)
            x1 = round(opt.centre(1)-w/2);
            y1 = round(opt.centre(2)-h/2);
            x2 = round(opt.centre(1)+w/2);
            y2 = round(opt.centre(2)+h/2);
        elseif ~isempty(opt.topleft)
            x1 = opt.topleft(1);
            y1 = opt.topleft(2);
            x2 = x1 + w;
            y2 = x1 + h;
        else
            error('must specify top left or centre');
        end
    else
       if ~isempty(opt.matlab)
            if numrows(opt.matlab) > 1
                for i=1:numrows(opt.matlab)
                    plot_box('matlab', opt.matlab(i,:), args{:});
                end
                return
            else
            x1 = opt.matlab(1);
            y1 = opt.matlab(2);
            x2 = opt.matlab(1) + opt.matlab(3);
            y2 = opt.matlab(2) + opt.matlab(4);
            end
       elseif all(size(args{1}) == [2 2])
            % first arg is a box
            b = args{1};
            x1 = b(1); y1 = b(2);
            x2 = b(3); y2 = b(4);
            args = args(2:end);
        else
            % use first 4 args as x1 y1 x2 y2
            x1 = args{1};
            y1 = args{2};
            x2 = args{3};
            y2 = args{4};
            args = args(5:end);
        end
    end
    x = [x1 x2 x2 x1 x1];
    y = [y1 y1 y2 y2 y1];
    
    
    holdon = ishold;
    hold on
    
    if isempty(opt.fillcolor)
        % outline only
        hh = plot(x, y, args{:})
    else
        % filled shape
        hh = patch(x, y, 0*y, 'FaceColor', opt.fillcolor, ...
            'FaceAlpha', opt.alpha, 'EdgeColor', opt.edgecolor, args{:});
    end

    if nargout > 0
        h = hh;
    end
    
    if holdon == 0
        hold off
    end
