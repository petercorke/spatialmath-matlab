%CIRCLE Compute points on a circle
%
% CIRCLE(C, R, OPTIONS) plots a circle centred at C (1x2) with radius R on the current
% axes.
%
% X = CIRCLE(C, R, OPTIONS) is a matrix (2xN) whose columns define the 
% coordinates [x,y] of points around the circumference of a circle 
% centred at C (1x2) and of radius R.
%
% C is normally 2x1 but if 3x1 then the circle is embedded in 3D, and X is Nx3.
% The circle is always in the xy-plane with a z-coordinate of C(3).
%
% Options::
%  'n',N   Specify the number of points (default 50)

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
