%PLOT_ELLIPSE Draw an ellipse on the current plot
%
% PLOT_ELLIPSE(A, LS) draws an ellipse defined by X'AX = 0 on the 
% current plot, centred at the origin, with Matlab line style LS.
%
% PLOT_ELLIPSE(A, C, LS) as above but centred at C=[X,Y].
% current plot.  If C=[X,Y,Z] the ellipse is parallel to the XY plane
% but at height Z.
%
% See also PLOT_CIRCLE.

function h = plot_ellipse(A, xc, varargin)

    if size(A,1) ~= size(A,2)
        error('ellipse is defined by a square matrix');
    end

    if size(A,1) > 3
        error('can only plot ellipsoid for 2 or 3 dimenions');
    end
    
    if nargin < 2
        xc = zeros(1, size(A,1));
    end
    if nargin < 3
        varargin = {};
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
            pe = bsxfun(@plus, xc(:), pe);
        end

        % put back to mesh format
        Xe = reshape(pe(1,:), size(Xs));
        Ye = reshape(pe(2,:), size(Ys));
        Ze = reshape(pe(3,:), size(Zs));

        % plot it
        h = mesh(Xe, Ye, Ze);

    else
        %% plot an ellipse
        
        opt.z = 0;
        [opt,varargin] = tb_optparse(opt, varargin);

        [V,D] = eig(A);
        
        % define points on a unit circle
        th = linspace(0, 2*pi, 50);
        pc = [cos(th);sin(th)];
        
        % warp it into the ellipse
        pe = sqrtm(A)*pc;

        % offset it to optional non-zero centre point
        if nargin > 1
            pe = bsxfun(@plus, xc(:), pe);
        end

        % plot it
        if opt.z == 0
            h = plot(pe(1,:), pe(2,:), varargin{:} );
        else
            h = plot3(pe(1,:), pe(2,:), ones(1,numcols(pe))*opt.z, varargin{:} );
            
        end
    end
    holdon = ishold;
    hold on
end
