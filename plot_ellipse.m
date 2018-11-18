%PLOT_ELLIPSE Draw an ellipse or ellipsoid
%
% PLOT_ELLIPSE(E, OPTIONS) draws an ellipse or ellipsoid defined by X'EX =
% 0 on the current plot, centred at the origin.  E (2x2) for an ellipse and
% E (2x3) for an ellipsoid.
%
% PLOT_ELLIPSE(E, C, OPTIONS) as above but centred at C=[X,Y].  If
% C=[X,Y,Z] the ellipse is parallel to the XY plane but at height Z.
%
% H = PLOT_ELLIPSE(E, C, OPTIONS) as above but return graphic handle.
%
% Animation::
%
% First draw the ellipse and keep its graphic handle, then alter it, eg.
%
%          H = PLOT_ELLIPSE(E, C, 'r')
%          PLOT_ELLIPSE(C, R, 'alter', H);
%
% Options::
% 'confidence',C   confidence interval, range 0 to 1
% 'alter',H        alter existing ellipses with handle H
% 'npoints',N      use N points to define the ellipse (default 40)
% 'edgecolor'      color of the ellipse boundary edge, MATLAB color spec
% 'fillcolor'      the color of the circle's interior, MATLAB color spec
% 'alpha'          transparency of the fillcolored circle: 0=transparent, 1=solid
% 'shadow'         show shadows on the 3 walls of the plot box
%
% - For an unfilled ellipse any standard MATLAB LineStyle such as 'r' or 'b---'.
% - For an unfilled ellipse any MATLAB LineProperty options can be given such as 'LineWidth', 2.
% - For a filled ellipse any MATLAB PatchProperty options can be given.
%
% Notes::
% - If A (2x2) draw an ellipse, else if A(3x3) draw an ellipsoid.
% - The ellipse is added to the current plot irrespective of hold status.
% - Shadow option only valid for ellipsoids.
% - If a confidence interval is given the scaling factor is com;uted using
%   an approximate inverse chi-squared function.
%
% See also PLOT_ELLIPSE_INV, PLOT_CIRCLE, PLOT_BOX, PLOT_POLY.


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

function handles = plot_ellipse(E, varargin)
    
    assert(size(E,1) == size(E,2), 'ellipse is defined by a square matrix');
    assert( size(E,1) == 2 || size(E,1) == 3, 'can only plot ellipsoid for 2 or 3 dimenions');
    
    opt.fillcolor = 'none';
    opt.alpha = 1;
    opt.edgecolor = 'k';
    opt.alter = [];
    opt.npoints = 40;
    opt.shadow = false;
    opt.confidence = [];
    
    [opt,arglist,ls] = tb_optparse(opt, varargin);

    % process some arguments
    
    if ~isempty(ls)
        opt.edgecolor = ls{1};
    end
    
    % process the probability
    if isempty(opt.confidence)
        s = 1;
    else 
        s = sqrt(chi2inv_rtb(opt.confidence, 2));
    end
    
    if length(arglist) > 0 && isnumeric(arglist{1})
        % ellipse centre is provided
        centre = arglist{1};
        arglist = arglist(2:end);
    else
        % default to origin
        centre = zeros(1, size(E,1));
    end

    % check the ellipse to be altered
    if ~isempty(opt.alter) & ~ishandle(opt.alter)
        error('RTB:plot_circle:badarg', 'argument to alter must be a valid graphic object handle');
    end
    
    holdon = ishold();
    hold on
    
    if size(E,1) == 3
        %% plot an ellipsoid
        
        % define mesh points on the surface of a unit sphere
        [Xs,Ys,Zs] = sphere();
        ps = [Xs(:) Ys(:) Zs(:)]';
        
        % warp it into the ellipsoid
        pe = sqrtm(E) * ps;
        
        % offset it to optional non-zero centre point
        if nargin > 1
            pe = bsxfun(@plus, centre(:), pe);
        end
        
        % put back to mesh format
        Xe = reshape(pe(1,:), size(Xs));
        Ye = reshape(pe(2,:), size(Ys));
        Ze = reshape(pe(3,:), size(Zs));
        

        if isempty(opt.alter)
              % plot it
%             Ce = ones(size(Xe));
%             Ce = cat(3, Ce*0.8, Ce*0.4, Ce*0.4);
            h = mesh(Xe, Ye, Ze, 'FaceColor', opt.fillcolor, ...
                        'FaceAlpha', opt.alpha, 'EdgeColor', opt.edgecolor, arglist{:});
        else
            % update an existing plot
            set(opt.alter, 'xdata', Xe, 'ydata', Ye, 'zdata', Ze,  ...
                        arglist{:});
        end
        
        % draw the shadow
        if opt.shadow
            I = ones(size(Xe));
            a = [xlim ylim zlim];
            mesh(a(1)*I, Ye, Ze, 'FaceColor', 0.7*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            mesh(Xe, a(3)*I, Ze, 'FaceColor', 0.7*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            mesh(Xe, Ye, a(5)*I, 'FaceColor', 0.7*[1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
        
    else
        %% plot an ellipse
        
                
        [V,D] = eig(E);
        
        % define points on a unit circle
        th = linspace(0, 2*pi, opt.npoints);
        pc = [cos(th);sin(th)];
        
        % warp it into the ellipse
        pe = sqrtm(E)*pc * s;
        
        % offset it to optional non-zero centre point
        centre = centre(:);
        if nargin > 1
            pe = bsxfun(@plus, centre(1:2), pe);
        end
        x = pe(1,:); y = pe(2,:);

%         if length(centre) > 2
%             % plot 3D data
%             z = ones(size(x))*centre(3);
%             if isempty(opt.alter)
%                 h = plot3(x', y', z', varargin{:});
%             else
%                 set(opt.alter, 'xdata', x, 'ydata', y, 'zdata', z, arglist{:});
%             end

            % plot 2D data

            if length(centre) > 2
                % plot 3D data
                z = ones(size(x))*centre(3);
            else
                z = zeros(size(x));
            end
            
            
            if strcmpi(opt.fillcolor, 'none')
                % outline only, draw a line
                
                if isempty(ls)
                    if ~isempty(opt.edgecolor)
                        arglist = ['Color', opt.edgecolor, arglist];
                    end
                else
                    arglist = [ls arglist];
                end

                if isempty(opt.alter)
                    h = plot3(x', y', z', arglist{:});
                else
                    set(opt.alter, 'xdata', x, 'ydata', y);
                end
            else
                % fillcolored, use a patch
                
                if ~isempty(opt.edgecolor)
                    arglist = ['EdgeColor', opt.edgecolor, arglist];
                end
                
                arglist = [ls, 'FaceAlpha', opt.alpha, arglist];
                
                                
                if isempty(opt.alter)
                    h = patch(x', y', z', opt.fillcolor, arglist{:});
                else
                    set(opt.alter, 'xdata', x, 'ydata', y);
                end
                
            end
        end
    
  if ~holdon
      hold off
  end
    
    if nargout > 0
        handles = h;
    end
end
