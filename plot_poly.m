%PLOT_POLY Plot a polygon
%
% PLOTPOLY(P, OPTIONS) plot a polygon defined by columns of P which
% can be 2xN or 3xN.
%
% OPTIONS::
%  'fill'    the color of the circle's interior, Matlab color spec
%  'alpha'   transparency of the filled circle: 0=transparent, 1=solid.
%
% See also PLOT, PATCH, Polygon.

% TODO: options for fill, not filled, line style, labels (cell array of strings)

function h_ = plot_poly(p, varargin)

    if numcols(p) < 3,
        error('too few points for a polygon');
    end

    opt.fill = [];
    opt.alpha = 1;

    [opt,arglist] = tb_optparse(opt, varargin);

    % default marker style
    if isempty(arglist)
        arglist = {'r-'};
    end

    ish = ishold();
	hold on

    x = [p(1,:) p(1,1)];
    y = [p(2,:) p(2,1)];
    if numrows(p) == 2
        % plot 2D data
        h(1) = plot(x, y, arglist{:});
        if ~isempty(opt.fill)
            h(2) = patch(x', y', 0*y', 'FaceColor', opt.fill, ...
                'FaceAlpha', opt.alpha);
        end
    elseif numrows(p) == 3
        % plot 3D data
        z = [p(3,:) p(3,1)];
        h(1) = plot3(x, y, z, arglist{:});
        if ~isempty(opt.fill)
            h(2) = patch(x, y, z, 0*y, 'FaceColor', opt.fill, ...
                'FaceAlpha', opt.alpha);
        end
    else
        error('point data must have 2 or 3 rows');
    end

    if ~ish
        hold off
    end
    %figure(gcf)
    
    if nargout > 0
        h_ = h;
    end
