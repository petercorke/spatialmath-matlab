%PLOT_POLY Draw a polygon
%
% plot_poly(P, OPTIONS) adds a closed polygon defined by vertices in the columns
% of P (2xN), in the current plot.
%
% H = plot_poly(...) as above but returns a graphics handle.
%
% plot_poly(H, )
%
% OPTIONS::
%  'fillcolor',F    the color of the circle's interior, MATLAB color spec
%  'alpha',A        transparency of the filled circle: 0=transparent, 1=solid.
%  'edgecolor',E    edge color
%  'animate'        the polygon can be animated
%  'tag',T          the polygon is created with a handle graphics tag
%  'axis',h         handle of axis or UIAxis to draw into (default is current axis)
%
% - For an unfilled polygon:
%   - any standard MATLAB LineStyle such as 'r' or 'b---'.
%   - any MATLAB LineProperty options can be given such as 'LineWidth', 2.
% - For a filled polygon any MATLAB PatchProperty options can be given.
%
% Notes::
% - If P (3xN) the polygon is drawn in 3D
% - If not filled the polygon is a line segment, otherwise it is a patch
%   object.
% - The 'animate' option creates an hgtransform object as a parent of the
%   polygon, which can be animated by the last call signature above.
% - The graphics are added to the current plot.
%
% Example::
%
%          POLY = [0 1 2; 0 1 0];
%          H = plot_poly(POLY, 'animate', 'r'); % draw a red polygon
%
%          H = plot_poly(POLY, 'animate', 'r'); % draw a red polygon that can be animated
%          plot_poly(H, transl(2,1,0) );  % transform its vertices by (2,1)
%
% See also PLOT_BOX, PLOT_CIRCLE, PATCH, Polygon.

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

% TODO: options for fill, not filled, line style, labels (cell array of strings)
%  handle filled/unfilled better, 'none' is synonymous with []?
%  is moveable used anywhere, seems broken
% TODO: move this animate logic to circle + ellipse
% TODO: move this 'axis' logic to circle + ellipse


function h_ = plot_poly(p, varargin)

    if ishandle(p)
        % PLOT_POLY(H, T)
        %  - animate existing polygon
        tr = varargin{1};
        if isvec(tr)
            m = SE2(tr); m = m.SE3; m = m.double;
        elseif ishomog2(tr)
            m = SE2(tr); m = m.SE3; m = m.double;
        elseif ishomog(tr)
            m = tr;
        else
            error('unknown transform type');
        end
        % set the transformation for the handle
        set(p, 'Matrix', m);
        return
    end
    
    % create a new polygon
    

    % process options
    opt.fillcolor = [];
    opt.fill = '@fillcolor';
    opt.alpha = 1;
    opt.animate = false;
    opt.edgecolor = 'None';
    opt.tag = [];
    opt.axis = [];
    
    [opt,args,ls] = tb_optparse(opt, varargin);


    if ~isempty(opt.fillcolor) && strcmp(opt.edgecolor, 'None')
        opt.edgecolor = 'k';
    end
    if isempty(opt.axis)
        opt.axis = gca;
    end
    

    % unpack the data and wrap it around to form a closed polygon
    assert( numcols(p) > 2, 'too few points for a polygon');
    assert( numrows(p) == 2 || numrows(p) == 3, 'data must have 2 or 3 rows');
    
    x = p(1,:); y = p(2,:);
    if numrows(p) == 3
        z = p(3,:);
    end
    
    if isempty(opt.fillcolor)
        % wrap it around to form a closed polygon
        x = [x x(1)];
        y = [y y(1)];
        if numrows(p) == 3
            z = [p(3,:) p(3,1)];
        end
    end
     
    if opt.animate
        if ~isempty(opt.tag)
            hg = hgtransform(opt.axis, 'Tag', opt.tag);
        else
            hg = hgtransform(opt.axis);
        end
        args = [args, 'Parent', {hg}];
    end

    if ~isempty(opt.fillcolor) && ~isempty(opt.edgecolor)
        % in fill mode, optionally set edge color
        args = [args, 'EdgeColor', opt.edgecolor];
    end

    ish = ishold();
	hold on

    switch numrows(p)
        case 2
            % plot 2D data
            if isempty(opt.fillcolor)
                h = plot(x, y, ls{:}, args{:});
            else
                h = patch(x', y', opt.fillcolor, ...
                    'EdgeColor', opt.edgecolor, ...
                    'FaceAlpha', opt.alpha, args{:});
            end
            
        case 3
            % plot 3D data
            if isempty(opt.fillcolor)
                h = plot3(x, y, z, args{:});
            else
                h = patch(x, y, z, opt.fillcolor, ...
                    'FaceAlpha', opt.alpha, args{:});
            end
    end

    if ~ish
        hold off
    end
    %figure(gcf)
    
    if nargout > 0
        if opt.animate
            h_ = hg;
        else
            h_ = h;
        end
    end
    
    
    %     if ~isempty(opt.handle)
%         if numrows(p) == 2
%             set(opt.handle, 'Xdata', x, 'Ydata', y);
%         elseif numrows(p) == 3
%             set(opt.handle, 'Xdata', x, 'Ydata', y, 'Zdata', z);
%         end
%         return;
%     end
    

    % default marker style
% %     if isempty(arglist)
% %         arglist = {'r-'};
% %     end
