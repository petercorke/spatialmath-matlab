%TRINTERP2 Interpolate SE(2) homogeneous transformations
%
% R = TRINTERP2(R0, R1, S) is a rotation matrix (2x2) interpolated
% between R0 when S=0 and R1 when S=1.  R0 and R1 are both rotation matrices
% (2x2).  If S (Nx1) then T (2x2xN) is a sequence of
% rotation matrices corresponding to the interpolation values in S.
%
% R = TRINTERP2(R1, S) as above but interpolated between the identity matrix
% when S=0 to T1 when S=1.
%
% T = TRINTERP2(T0, T1, S) is a homogeneous transform (3x3) interpolated
% between T0 when S=0 and T1 when S=1.  T0 and T1 are both homogeneous
% transforms (3x3).  If S (Nx1) then T (3x3xN) is a sequence of
% homogeneous transforms corresponding to the interpolation values in S.
%
% T = TRINTERP2(T1, S) as above but interpolated between the identity matrix
% when S=0 to T1 when S=1.
%
% See also trinterp, SE3.interp, UnitQuaternion.




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

function T = trinterp2(A, B, C)

    switch nargin
        case 2
            % trinterp(T, s)
            T1 = A; s = B;
            
            th0 = 0;
            th1 = atan2(T1(2,1), T1(1,1));
            if ~isrot2(T1)
                p0 = [0 0]';
                p1 = transl2(T1);
            end
        case 3
            % trinterp(T1, T2, s)
            T0 = A; T1 = B; s = C;
            assert(all(size(A) == size(B)), 'RTB:trinterp2:badarg', '2 matrices must be same size');
            th0 = atan2(T0(2,1), T0(1,1));
            th1 = atan2(T1(2,1), T1(1,1));
            if ~isrot2(T0)
                p0 = transl2(T0);
                p1 =transl2(T1);
            end
        otherwise
            error('RTB:trinterp2:badarg', 'must be 2 or 3 arguments');
    end
    
    if length(s) == 1 && s > 1 && (s == floor(s))
        % integer value
        s = [0:(s-1)] / (s-1);
    elseif any(s<0 | s>1)
        error('RTB:trinterp2:badarg', 'values of S outside interval [0,1]');
    end
    
    if isrot2(T1)
        
        % SO(2) case
        for i=1:length(s)
            th = th0*(1-s(i)) + s(i)*th1;
            
            T(:,:,i) = rot2(th);
        end
    else
        % SE(2) case
        for i=1:length(s)
            th = th0*(1-s(i)) + s(i)*th1;
            pr = p0*(1-s(i)) + s(i)*p1;
            
            T(:,:,i) = rt2tr(rot2(th), pr);
        end
    end
    
end
