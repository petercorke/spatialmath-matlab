%RT2TR Convert rotation and translation to homogeneous transform
%
% TR = RT2TR(R, t) is a homogeneous transformation matrix (N+1xN+1) formed
% from an orthonormal rotation matrix R (NxN) and a translation vector t
% (Nx1).  Works for R in SO(2) or SO(3):
%  - If R is 2x2 and t is 2x1, then TR is 3x3
%  - If R is 3x3 and t is 3x1, then TR is 4x4
%
% For a sequence R (NxNxK) and t (NxK) results in a transform sequence (N+1xN+1xK).
%
% Notes::
% - The validity of R is not checked
%
% See also T2R, R2T, TR2RT.



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

function T = rt2tr(R, t)

    assert( size(R,1) == size(R,2), 'RTB:rt2tr:badarg', 'R must be square');
    
    n = size(R,2);
    B = [zeros(1,n) 1];
    
    if size(R,3) > 1
        % vector case
        assert( size(R,1) == size(t,1), 'RTB:rt2tr:badarg', 'R and t must have conforming dimensions')
        assert( size(R,3) == size(t,2), 'RTB:rt2tr:badarg', 'For sequence size(R,3) must equal size(t,2)');
        
        T = zeros(n+1,n+1,size(R,3));
        for i=1:size(R,3)
            T(:,:,i) = [R(:,:,i) t(:,i); B];
        end
    else
        % scalar case
        assert( isvec(t, size(R,1)), 'RTB:rt2tr:badarg', 'R and t must have conforming dimensions')
        T = [R t(:); B];
    end
end

