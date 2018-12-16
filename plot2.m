%PLOT2 Plot trajectories
%
% Convenience function for plotting 2D or 3D trajectory data stored in a
% matrix with one row per time step.
%
% PLOT2(P) plots a line with coordinates taken from successive rows of P.  P
% can be Nx2 or Nx3.
%
% If P has three dimensions, ie. Nx2xM or Nx3xM then the M trajectories are
% overlaid in the one plot.
%
% PLOT2(P, LS) as above but the line style arguments LS are passed to plot.
%
% See also MPLOT, PLOT.


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
function h = plot2(p1, varargin)

    if ndims(p1) == 2
        switch numcols(p1)
            case 3
                hh = plot3(p1(:,1), p1(:,2), p1(:,3), varargin{:});
            case 2
                hh = plot(p1(:,1), p1(:,2), varargin{:});
            otherwise
                error('RTB:plot2:badarg', 'Data must have 2 or 3 columns');
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
