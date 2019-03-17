%PLOT2 Plot trajectories
%
% Convenience function for plotting 2D or 3D trajectory data stored in a
% matrix with one row per time step.
%
% PLOT2(P) plots a line with coordinates taken from successive rows of P(Nx2 or Nx3).
%
% If P has three dimensions, ie. Nx2xM or Nx3xM then the M trajectories are
% overlaid in the one plot.
%
% PLOT2(P, LS) as above but the line style arguments LS are passed to plot.
%
% See also MPLOT, PLOT.

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
function h = plot2(p1, varargin)

    if ndims(p1) == 2
        switch numcols(p1)
            case 3
                hh = plot3(p1(:,1), p1(:,2), p1(:,3), varargin{:});
            case 2
                hh = plot(p1(:,1), p1(:,2), varargin{:});
            otherwise
                error('SMTB:plot2:badarg', 'Data must have 2 or 3 columns');
        end
        if nargout == 1
            h = hh;
        end
    else
        clf
        hold on
        for i=1:size(p1,2)
            plot2( squeeze(p1(:,:,i)) );
        end
    end
