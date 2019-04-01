%PLOT_RIBBON Draw a wide curved 3D arrow
%
% plot_ribbon() adds a 3D curved arrow "ribbon" to the current plot.  The ribbon by
% default is about the z-axis at the origin.
%
% Options::
% 'radius',R     radius of the ribbon (default 0.25)
% 'N',N          number of points along the ribbon (default 100)
%
% 'd',D          ratio of shaft length to total  (default 0.9)
% 'w1',W         width of shaft (default 0.2)
% 'w2',W         width of head (default 0.4)
% 'phi',P        length of ribbon as fraction of circle (default 0.8)
% 'phase',P      rotate the arrow about its axis (radians, default 0)
%
% 'color',C      color as MATLAB ColorSpec (default 'r')
% 'specular',S   specularity of surface (default 0.2)
% 'diffuse',D    diffusivity of surface (default 0.8)
%
% 'nice'         adjust the phase for nicely phased arrow 
%
% The parameters of the ribbon are:
%
%          ^
%          |                            | \
%          |  ^  +----------------------|  \
%          |  |  |                          .
%          |  v  +----------------------|  /
%          |  w1                        | /
%          v     <---------- d --------->
%         w2    <----------- phi ---------->
%
% Examples::
%
% To draw the ribbon at distance A along the X, Y, Z axes is:
%          plot_ribbon2( SE3(A,0,0)*SE3.Ry(pi/2) )
%          plot_ribbon2( SE3(0, A,0)*SE3.Rx(pi/2) )
%          plot_ribbon2( SE3(0, 0, A) )
%          shading interp
%          camlight
%
% See also plot_arrow, plot.


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

% clf
% trplot( SE3 )


function plot_ribbon2(pose, varargin)

    opt.radius = 0.25;
    opt.N = 100;

    opt.d = 0.90;
    opt.w1 = 0.2;
    opt.w2 = 0.4;
    opt.phi = 0.8;
    opt.phase = 0;

    opt.nice = false;

    opt.color = 'r';
    opt.specular = 0.2;
    opt.diffuse = 0.8;

    opt = tb_optparse(opt, varargin);

    Nd = floor(opt.N*opt.d);

    if opt.nice
        % ensure that gap in the ribbon is towards the viewpoint
        c = circle([0 0 0], opt.radius, 'n', opt.N);

        p = get(gca, 'CameraPosition')';

        [~,k] = min( colnorm( pose*c-p ) );
        opt.phase = opt.phase + (k - 10) / opt.N * 2*pi;
    end
    
    % compute canonic arrow about the z-axis and centered at origin
    
    theta = linspace(0, opt.phi*2*pi, opt.N) + opt.phase + opt.phi*2*pi/2;
    % replicate the Nd'th point, this means that
    % x(Nd) = x(Nd+1)
    % y(Nd) = y(Nd+1)
    % which makes the base of the arrowhead normal to the shaft
    theta = [theta(1:Nd) theta(Nd) theta(Nd+1:end)];
    x = opt.radius * cos(theta);
    y = opt.radius * sin(theta);
    
    % compute the width of the arrow, varies along length
    T = [0 opt.w1/2; Nd-1 opt.w1/2; Nd opt.w2/2; opt.N-1 0];
    z = interp1(T(:,1), T(:,2), [0:opt.N], 'linear');
    
    
    % build the mesh matrices
    X = [x; x];
    Y = [y; y];
    Z = [-z; z];
    
    z1 = ones(size(Z));
    C = cat(3, z1*1, z1*0, z1*0);
    hg = hgtransform;

    ish = ishold()
    hold on
    surf(X, Y, Z, C, 'Parent', hg, 'EdgeColor', 'None', 'FaceColor', opt.color, ...
        'BackFaceLighting', 'reverselit', 'SpecularStrength', opt.specular, 'DiffuseStrength', opt.diffuse);
    hold (ish)
    
    hg.Matrix = pose.T;
end