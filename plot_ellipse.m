%PLOT_ELLIPSE Draw an ellipse on the current plot
%
% PLOT_ELLIPSE(A, LS) draws an ellipse defined by X'AX = 0 on the
% current plot, centred at the origin, with Matlab line style LS.
%
% PLOT_ELLIPSE(A, C, LS) as above but centred at C=[X,Y].
% current plot.  If C=[X,Y,Z] the ellipse is parallel to the XY plane
% but at height Z.
%
% H = PLOT_CIRCLE(C, R, options) as above but return handles. For multiple
% circles H is a vector of handles, one per circle.
%
% Options::
% 'edgecolor'   the color of the circle's edge, Matlab color spec
% 'fillcolor'   the color of the circle's interior, Matlab color spec
% 'alpha'       transparency of the filled circle: 0=transparent, 1=solid
% 'alter',H     alter existing circles with handle H
%
% See also PLOT_CIRCLE.

function handles = plot_ellipse(A, centre, varargin)
    
    if size(A,1) ~= size(A,2)
        error('ellipse is defined by a square matrix');
    end
    
    if size(A,1) > 3
        error('can only plot ellipsoid for 2 or 3 dimenions');
    end
    
    if nargin < 2
        centre = zeros(1, size(A,1));
    end
    if nargin < 3
        varargin = {};
    end
    
    
    opt.fillcolor = [];
    opt.alpha = 1;
    opt.edgecolor = 'k';
    opt.alter = [];
    
    [opt,arglist] = tb_optparse(opt, varargin);
    
    if ~isempty(opt.alter) & ~ishandle(opt.alter)
        error('RTB:plot_circle:badarg', 'argument to alter must be a valid graphic object handle');
    end
    
    holdon = ishold();
    hold on
    
    if size(A,1) == 3
        %% plot an ellipsoid
        
        % define mesh points on the surface of a unit sphere
        [Xs,Ys,Zs] = sphere();
        ps = [Xs(:) Ys(:) Zs(:)]';
        
        % warp it into the ellipsoid
        pe = sqrtm(A) * ps;
        
        % offset it to optional non-zero centre point
        if nargin > 1
            pe = bsxfun(@plus, centre(:), pe);
        end
        
        % put back to mesh format
        Xe = reshape(pe(1,:), size(Xs));
        Ye = reshape(pe(2,:), size(Ys));
        Ze = reshape(pe(3,:), size(Zs));
        
        % plot it
        if isempty(opt.alter)
            h = mesh(Xe, Ye, Ze, arglist{:});
        else
            set(opt.alter, 'xdata', Xe, 'ydata', Ye, 'zdata', Ze, arglist{:});
            
        end
        
    else
        %% plot an ellipse
        
        
        [V,D] = eig(A);
        
        % define points on a unit circle
        th = linspace(0, 2*pi, 50);
        pc = [cos(th);sin(th)];
        
        % warp it into the ellipse
        pe = sqrtm(A)*pc;
        
        % offset it to optional non-zero centre point
        centre = centre(:);
        if nargin > 1
            pe = bsxfun(@plus, centre(1:2), pe);
        end
        x = pe(1,:); y = pe(2,:);
        

        if length(centre) > 2
            % plot 3D data
            z = ones(size(x))*centre(3);
            if isempty(opt.alter)
                h = plot3(x, y, z, varargin{:});
            else
                set(opt.alter, 'xdata', x, 'ydata', y, 'zdata', z, arglist{:});
            end
        else
            % plot 2D data
            if isempty(opt.fillcolor)
                if isempty(opt.alter)
                    h = plot(x, y, arglist{:});
                else
                    set(opt.alter, 'xdata', x, 'ydata', y, arglist{:});
                end
            else
                if isempty(opt.alter)
                    h = patch(x, y, 0*y, 'FaceColor', opt.fillcolor, ...
                        'FaceAlpha', opt.alpha, 'EdgeColor', opt.edgecolor, arglist{:});
                else
                    set(opt.alter, 'xdata', x, 'ydata', y, arglist{:});
                end
                
            end
        end
    end
    holdon = ishold;
    hold on
    
    if nargout > 0
        handles = h;
    end
end
