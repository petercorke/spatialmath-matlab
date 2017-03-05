%PLOT_POLY Draw a polygon
%
% PLOT_POLY(P, OPTIONS) adds a polygon defined by columns of P (2xN), in the current
% plot with default line style.
%
% H = PLOT_POLY(P, OPTIONS) as above but processes additional options and
% returns a graphics handle.
%
% Animation::
%
% PLOT_POLY(H, T) sets the pose of the polygon with handle H to the pose
% given by T (3x3 or 4x4).
%
% Create a polygon that can be animated, then alter it, eg.
%
%          H = PLOT_POLY(P, 'animate', 'r')
%          PLOT_POLY(H, transl(2,1,0) );
%
% OPTIONS::
%  'fillcolor',F  the color of the circle's interior, MATLAB color spec
%  'alpha',A      transparency of the filled circle: 0=transparent, 1=solid.
%  'edgecolor',E  edge color
%  'animate'      the polygon can be animated
%  'tag',T        the polygon is created with a handle graphics tag
%
% - For an unfilled polygon any standard MATLAB LineStyle such as 'r' or 'b---'.
% - For an unfilled polygon any MATLAB LineProperty options can be given such as 'LineWidth', 2.
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
% See also PLOT_BOX, PLOT_CIRCLE, PATCH, Polygon.


% Copyright (C) 1993-2017, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for MATLAB (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.
%
% http://www.petercorke.com

% TODO: options for fill, not filled, line style, labels (cell array of strings)
%  handle filled/unfilled better, 'none' is synonymous with []?
%  is moveable used anywhere, seems broken
% TODO: move this animate logic to circle + ellipse


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
    opt.fill = [];
    opt.alpha = 1;
    opt.animate = false;
    opt.edgecolor = 'None';
    opt.tag = [];
    
    [opt,args,ls] = tb_optparse(opt, varargin);

        
    if ~isempty(opt.fill)
        opt.fillcolor = opt.fill;
    end
    if ~isempty(opt.fillcolor) && strcmp(opt.edgecolor, 'None')
        opt.edgecolor = 'k';
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
            hg = hgtransform('Tag', opt.tag);
        else
            hg = hgtransform();
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
