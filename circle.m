%CIRCLE Compute points on a circle
%
% CIRCLE(C, R, OPT) plot a circle centred at C with radius R.
%
% X = CIRCLE(C, R, OPT) return an Nx2 matrix whose rows define the 
% coordinates [x,y] of points around the circumferance of a circle 
% centred at C and of radius R.
%
% C is normally 2x1 but if 3x1 then the circle is embedded in 3D, and X is Nx3,
% but the circle is always in the xy-plane with a z-coordinate of C(3).
%
% Options::
%  'n',N   Specify the number of points (default 50)
function out = circle(centre, rad, varargin)

	opt.n = 50;
    
    [opt,arglist] = tb_optparse(opt, varargin);

    % compute points on circumference
	th = [0:opt.n-1]'/ opt.n*2*pi;
    x = rad*cos(th) + centre(1);
    y = rad*sin(th) + centre(2);

    % add extra row if z-coordinate is specified, but circle is always in xy plane
    if length(centre) > 2
        z = ones(size(x))*centre(3);
        p = [x y z]';
    else
        p = [x y]';
    end

    if nargout > 0
        % return now
        out = p;
        return;
    else
        % else plot the circle
        p = [p p(:,1)]; % make it close
        if numrows(p) > 2
            plot3(p(1,:), p(2,:), p(3,:), arglist{:});
        else
            plot(p(1,:), p(2,:), arglist{:});
        end
    end
