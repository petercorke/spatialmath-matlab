%PLOT_POLY Draw a polygon
%
% PLOT_POLY(P, OPTIONS) draws a polygon defined by columns of P (2xN), in the current plot.
%
% OPTIONS::
%  'fill',F    the color of the circle's interior, MATLAB color spec
%  'alpha',A   transparency of the filled circle: 0=transparent, 1=solid.
%
% Notes::
% - If P (3xN) the polygon is drawn in 3D
% - The line(s) is added to the current plot.
%
% See also PLOT_BOX, PATCH, Polygon.

% Copyright (C) 1993-2014, by Peter I. Corke
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

function h_ = plot_poly(p, varargin)

    if ishandle(p)
        tr = varargin{1};
        if isvec(tr)
            m = SE2(tr); m = m.SE3; m = m.double;
        elseif ishomog(tr)
            m = tr;
        else
            error('unknown transform type');
        end
        set(p, 'Matrix', m);
        return
    end
    
    if numcols(p) < 3,
        error('too few points for a polygon');
    end
    
    % unpack the data and wrap it around to form a closed polygon
    x = [p(1,:) p(1,1)];
    y = [p(2,:) p(2,1)];
    if numrows(p) == 3
        z = [p(3,:) p(3,1)];
    end
     
    opt.fill = 'none';
    opt.alpha = 1;
    %opt.handle = [];  no longer supported
    opt.moveable = false;
    opt.edge = [];
    opt.tag = [];
    if isempty(opt.fill) && isempty(opt.edge)
        opt.edge = 'k';
    end

    [opt,arglist,ls] = tb_optparse(opt, varargin);
    
    if opt.moveable
        if ~isempty(opt.tag)
            hg = hgtransform('Tag', opt.tag);
        else
            hg = hgtransform();
        end
        arglist = [arglist, 'Parent', {hg}];
    end

    if ~isempty(opt.edge)
    arglist = [arglist, 'EdgeColor', opt.edge];
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

    ish = ishold();
	hold on


    if numrows(p) == 2
        % plot 2D data
        if strcmp(opt.fill, 'none')
            h = plot(x, y, ls{:}, arglist{:});
        else
            plot(x, y, ls{:}, arglist{:})
            h = patch(x', y', 0*y', 'FaceColor', opt.fill, ...
                'FaceAlpha', opt.alpha, arglist{:});
        end
    elseif numrows(p) == 3
        % plot 3D data
        if isempty(opt.fill)
                    h = plot3(x, y, z, arglist{:});

        else
            h = patch(x, y, z, 0*y, 'FaceColor', opt.fill, ...
                'FaceAlpha', opt.alpha, arglist{:});
        end
    else
        error('point data must have 2 or 3 rows');
    end

    if ~ish
        hold off
    end
    %figure(gcf)
    
    if nargout > 0
        if opt.moveable
            h_ = hg;
        else
            h_ = h;
        end
    end
