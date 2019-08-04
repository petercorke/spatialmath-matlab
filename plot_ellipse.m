%PLOT_ELLIPSE Draw an ellipse or ellipsoid
%
% plot_ellipse(E, OPTIONS) draws an ellipse or ellipsoid defined by X'EX =
% 0 on the current plot, centred at the origin.  E (2x2) for an ellipse and
% E (2x3) for an ellipsoid.
%
% plot_ellipse(E, C, OPTIONS) as above but centred at C=[X,Y].  If
% C=[X,Y,Z] the ellipse is parallel to the XY plane but at height Z.
%
% H = plot_ellipse(...) as above but return graphic handle.
%
% Options::
% 'confidence',C   confidence interval, range 0 to 1
% 'alter',H        alter existing ellipses with handle H
% 'npoints',N      use N points to define the ellipse (default 40)
% 'edgecolor'      color of the ellipse boundary edge, MATLAB color spec
% 'fillcolor'      the color of the ellipses's interior, MATLAB color spec
% 'alpha'          transparency of the fillcolored ellipse: 0=transparent, 1=solid
% 'shadow'         show shadows on the 3 walls of the plot box
%
% - For an unfilled ellipse:
%   - any standard MATLAB LineStyle such as 'r' or 'b---'.
%   - any MATLAB LineProperty options can be given such as 'LineWidth', 2.
% - For a filled ellipse any MATLAB PatchProperty options can be given.
%
% Example::
%
%          H = plot_ellipse(diag([1 2]), [3 4]', 'r'); % draw red ellipse
%          plot_ellipse(diag([1 2]), [5 6]', 'alter', H); % move the ellipse
%          plot_ellipse(diag([1 2]), [5 6]', 'alter', H, 'LineColor', 'k'); % change color
%
%          plot_ellipse(COVAR, 'confidence', 0.95); % draw 95% confidence ellipse
%
% Notes::
% - The 'alter' option can be used to create a smooth animation.
% - If E (2x2) draw an ellipse, else if E (3x3) draw an ellipsoid.
% - The ellipse is added to the current plot irrespective of hold status.
% - Shadow option only valid for ellipsoids.
% - If a confidence interval is given then E is interpretted as a covariance
%   matrix and the ellipse size is computed using an inverse chi-squared function.
%   This requires CHI2INV in the Statistics and Machine Learning Toolbox or
%   CHI2INV_RTB from the Robotics Toolbox for MATLAB.
%
% See also PLOT_ELLIPSE_INV, PLOT_CIRCLE, PLOT_BOX, PLOT_POLY, CH2INV.

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
        if exist('chi2inv') == 2
            s = sqrt(chi2inv(opt.confidence, 2));
        elseif exist('chi2inv_rtb') == 2
            s = sqrt(chi2inv_rtb(opt.confidence, 2));
        else
            error('SMTB:missingfunc', 'Requires Stats Toolbox or RTB to be installed');
        end
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
        error('SMTB:plot_ellipse:badarg', 'argument to alter must be a valid graphic object handle');
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
