%PLOT_SPHERE Draw sphere
%
% PLOT_SPHERE(C, R, LS) draws spheres in the current plot.  C is the 
% centre of the sphere (3x1), R is the radius and LS is an optional MATLAB 
% ColorSpec, either a letter or a 3-vector.  
%
% PLOT_SPHERE(C, R, COLOR, ALPHA) as above but ALPHA specifies the opacity
% of the sphere where 0 is transparant and 1 is opaque.  The default is 1.
%
% If C (3xN) then N sphhere are drawn and H is Nx1.  If R (1x1) then all
% spheres have the same radius or else R (1xN) to specify the radius of
% each sphere.
%
% H = PLOT_SPHERE(...) as above but returns the handle(s) for the
% spheres.
%
% Notes::
% - The sphere is always added, irrespective of figure hold state.
% - The number of vertices to draw the sphere is hardwired.
%
% Example::
%         plot_sphere( mkgrid(2, 1), .2, 'b'); % Create four spheres
%         lighting gouraud  % full lighting model
%         light
%
% See also: plot_point, plot_box, plot_circle, plot_ellipse, plot_poly.

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

% TODO
% inconsistant call format compared to other plot_xxx functions.

function out = plot_sphere(c, r, varargin)

    opt.color = 'b';
    opt.alpha = 1;
    opt.mesh = 'none';
    opt.n = 40;

    [opt,args] = tb_optparse(opt, varargin);
    
    % backward compatibility with RVC
    if length(args) > 0
        opt.color = args{1};
    end
    if length(args) > 1
        opt.alpha = args{2};
    end
    
    daspect([1 1 1])
    hold_on = ishold;
    hold on
    [xs,ys,zs] = sphere(opt.n);

    if isvec(c,3)
        c = c(:);
    end
    if size(r) == 1
        r = r * ones(numcols(c),1);
    end

    if nargin < 4
        alpha = 1;
    end

    % transform the sphere
    for i=1:numcols(c)
        x = r(i)*xs + c(1,i);
        y = r(i)*ys + c(2,i);
        z = r(i)*zs + c(3,i);
                
        % the following displays a nice smooth sphere with glint!
        h = surf(x,y,z, ones(size(z)), 'FaceColor', opt.color, 'EdgeColor', opt.mesh, 'FaceAlpha', opt.alpha);
        % camera patches disappear when shading interp is on
        %h = surfl(x,y,z)
    end
    %lighting gouraud
    %light
    if ~hold_on
        hold off
    end
    if nargout > 0
        out = h;
    end
