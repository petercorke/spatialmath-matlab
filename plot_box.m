%PLOT_BOX	Draw a box on the current plot
%
% PLOT_BOX(B, LS) draws a box defined by B=[XL XR; YL YR] with optional Matlab
% linestyle options LS.
%
% PLOT_BOX(X1,Y1, X2,Y2, LS) draws a box with corners at (X1,Y1) and (X2,Y2),
% and optional Matlab linestyle options LS.
%
% PLOT_BOX('centre', P, 'size', W, LS) draws a box with center at P=[X,Y] and
% with dimensions W=[WIDTH HEIGHT].
%
% PLOT_BOX('topleft', P, 'size', W, LS) draws a box with top-left at P=[X,Y] 
% and with dimensions W=[WIDTH HEIGHT].

% Options::
% 'size',SZ     Specify size of box.  SZ=[WIDTH, HEIGHT] or if scalar then
%               WIDTH=HEIGHT=SZ
% 'topleft',P   Specify position of box by top left coordinate P
% 'centre',P    Specify position of box by centre coordinate P
%
% Additional options LS are passed to PLOT.
%
% See also plot.

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

function plot_box(varargin)
    opt.centre = [];
    opt.topleft = [];
    opt.size = [];

    [opt,varargin] = tb_optparse(opt, varargin);

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
        if all(size(varargin{1}) == [2 2])
            % first arg is a box
            b = varargin{1};
            x1 = b(1); y1 = b(2);
            x2 = b(3); y2 = b(4);
            varargin = varargin(2:end);
        else
            % use first 4 args as x1 y1 x2 y2
            x1 = varargin{1};
            y1 = varargin{2};
            x2 = varargin{3};
            y2 = varargin{4};
            varargin = varargin(5:end);
        end
    end
    p = [	x1 y1
            x2 y1
            x2 y2
            x1 y2
            x1 y1 ];

    holdon = ishold;
    hold on

    plot(p(:,1), p(:,2), varargin{:})

    if holdon == 0
        hold off
    end
