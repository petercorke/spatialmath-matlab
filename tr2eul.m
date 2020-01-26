%TR2EUL Convert SO(3) or SE(3) matrix to Euler angles
%
% EUL = TR2EUL(T, OPTIONS) are the ZYZ Euler angles (1x3) corresponding to
% the rotational part of a homogeneous transform T (4x4). The 3 angles
% EUL=[PHI,THETA,PSI] correspond to sequential rotations about the Z, Y and
% Z axes respectively.
%
% EUL = TR2EUL(R, OPTIONS) as above but the input is an orthonormal
% rotation matrix R (3x3).
%
% If R (3x3xK) or T (4x4xK) represent a sequence then each row of EUL
% corresponds to a step of the sequence.
%
% Options::
%  'deg'      Compute angles in degrees (radians default)
%  'flip'     Choose first Euler angle to be in quadrant 2 or 3.
%
% Notes::
% - There is a singularity for the case where THETA=0 in which case PHI is arbitrarily
%   set to zero and PSI is the sum (PHI+PSI).
% - Translation component is ignored.
%
% See also  EUL2TR, TR2RPY.

%## 3d rotation homogeneous

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

function euler = tr2eul(RR, varargin)
    
    assert(isrot(RR) || ishomog(RR), 'SMTB:tr2eul:badarg', 'argument must be a 3x3 or 4x4 matrix');

    opt.deg = false;
    opt.flip = false;
    opt = tb_optparse(opt, varargin);
    
    euler = zeros(size(RR,3),3);
    
    for i=1:size(RR,3)
        
        R = RR(:,:,i);
        
        % Method as per Paul, p 69.
        % euler = [phi theta psi]
        %
        
        if abs(R(1,3)) < eps && abs(R(2,3)) < eps
            % singularity
            eul(1) = 0;
            sp = 0;
            cp = 1;
            eul(2) = atan2(cp*R(1,3) + sp*R(2,3), R(3,3));
            eul(3) = atan2(-sp * R(1,1) + cp * R(2,1), -sp*R(1,2) + cp*R(2,2));
        else
            % non singular
            
            % Only positive phi is returned.
            if opt.flip
                eul(1) = atan2(-R(2,3), -R(1,3));
            else
                eul(1) = atan2(R(2,3), R(1,3));
            end
            sp = sin(eul(1));
            cp = cos(eul(1));
            eul(2) = atan2(cp*R(1,3) + sp*R(2,3), R(3,3));
            eul(3) = atan2(-sp * R(1,1) + cp * R(2,1), -sp*R(1,2) + cp*R(2,2));
        end
        
        euler(i,:) = eul;
    end
    
    if opt.deg
        euler = euler * 180/pi;
    end
end
