%PLOT_ARROW Draw an arrow in 2D or 3D
%
% PLOT_ARROW(P1, P2, OPTIONS) draws an arrow from P1 to P2 (2x1 or 3x1).
%
% PLOT_ARROW(P, OPTIONS) as above where the columns of P (2x2 or 3x2) define where P=[P1 P2].
%
% Options::
% - All options are passed through to arrow3.  
% - MATLAB colorspec such as 'r' or 'b--'
%
% See also ARROW3.


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
function plot_arrow(p1, varargin)
    
    if min(size(p1)) == 1
        % p1 is a vector
        p1 = p1(:);
        p2 = varargin{1};
        p2 = p2(:);
        assert(numrows(p1) == numrows(p2), 'RTB:plot_arrow', 'P1 and P2 must be the same length');
        varargin = varargin{2:end};
    else
        % p1 is a 2-column matrix
        assert(numcols(p1) == 2, 'RTB:plot_arrow', 'P1 must have 2 columns');
        p2 = p1(:,2);
        p1 = p1(:,1);
    end
    
    assert(any(numrows(p1) == [2 3]), 'RTB:plot_arrow', '2D or 3D points only');

    arrow3(p1', p2', varargin{:});
