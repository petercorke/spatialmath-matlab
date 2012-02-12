%PLOT_CIRCLE Draw a circle on the current plot
%
% PLOT_CIRCLE(C, R, options) draws a circle on the current plot with 
% centre C=[X,Y] and radius R.  If C=[X,Y,Z] the circle is drawn in the
% XY-plane at height Z.
%
% Options::
%  'edgecolor'   the color of the circle's edge, Matlab color spec
%  'fillcolor'   the color of the circle's interior, Matlab color spec
%  'alpha'       transparency of the filled circle: 0=transparent, 1=solid.

% Copyright (C) 1995-2009, by Peter I. Corke
%
% This file is part of The Machine Vision Toolbox for Matlab (MVTB).
% 
% MVTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MVTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with MVTB.  If not, see <http://www.gnu.org/licenses/>.
function plot_circle(centre, rad, varargin)

    opt.fillcolor = [];
    opt.alpha = 1;
    opt.edgecolor = 'k';

    [opt,arglist] = tb_optparse(opt, varargin);

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
            plot3(x, y, z, varargin{:});
        else
            % plot 2D data
            if isempty(opt.fillcolor)
                plot(x, y, arglist{:});
            else
                patch(x, y, 0*y, 'FaceColor', opt.fillcolor, ...
                    'FaceAlpha', opt.alpha, 'EdgeColor', opt.edgecolor, arglist{:})
            end
        end
    end

    if holdon == 0
        hold off
    end
