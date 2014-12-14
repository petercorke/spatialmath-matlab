%plot_ellipse Plot an ellipse
%
%   plot_ellipse(A, xc, ls)
%
%   ls is the standard line styles.

% Copyright (C) 1993-2014, by Peter I. Corke
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

function h = plot_ellipse_inv(A, xc, varargin)

    if nargin == 1
        h = plot_ellipse(inv(A));
    elseif nargin == 2
        h = plot_ellipse(inv(A), xc);
    else
        h = plot_ellipse(inv(A), xc, varargin{:});
    end
