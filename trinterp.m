%TRINTERP Interpolate SE(3) homogeneous transformations
%
% TRINTERP(T0, T1, S) is a homogeneous transform (4x4) interpolated
% between T0 when S=0 and T1 when S=1.  T0 and T1 are both homogeneous
% transforms (4x4).  If S (Nx1) then T (4x4xN) is a sequence of
% homogeneous transforms corresponding to the interpolation values in S.
%
% TRINTERP(T1, S) as above but interpolated between the identity matrix
% when S=0 to T1 when S=1.
%
% TRINTERP(T0, T1, M) as above but M is a positive integer and return a
% sequence (4x4xM) of homogeneous transforms linearly interpolating between 
% T0 and T1 in M steps.
%
% TRINTERP(T1, M) as above but return a sequence (4x4xM) of
% homogeneous interpolating between identity matrix and T1 in M steps.
%
% Notes::
% - T0 or T1 can also be an SO(3) rotation matrix (3x3) in which case the
%   result is (3x3xN).
% - Rotation is interpolated using quaternion spherical linear interpolation (slerp).
% - To obtain smooth continuous motion S should also be smooth and continuous,
%   such as computed by tpoly or lspb. 
%
% See also TRINTERP2, CTRAJ, SE3.interp, UnitQuaternion, TPOLY, LSPB.

%## 3d homogeneous trajectory

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

function T = trinterp(A, B, C)
    
    if nargin == 3
        %	TRINTERP(T0, T1, s)
        T0 = A; T1 = B; s = C(:)';
        
        if length(s) == 1 && s > 1 && (s == floor(s))
            % TRINTERP(T0, T1, M)
            s = linspace(0, 1, s);
        end
        assert(all(s>=0 & s<=1), 'SMTB:trinterp:badarg', 'values of S outside interval [0,1]');
        
        q0 = UnitQuaternion(T0);
        q1 = UnitQuaternion(T1);
        
        p0 = transl(T0);
        p1 = transl(T1);
        
        for i=1:length(s)
            qr = q0.interp(q1, s(i));
            pr = p0*(1-s(i)) + s(i)*p1;
            T(:,:,i) = rt2tr(qr.R, pr);
        end
    elseif nargin == 2
        %	TRINTERP(T, s)
        T0 = A; s = B(:)';
        
        if length(s) == 1 && s > 1 && (s == floor(s))
            % TRINTERP(T0, T1, M)
            s = linspace(0, 1, s);
        elseif any(s<0 | s>1)
            error('SMTB:trinterp:badarg', 'values of S outside interval [0,1]');
        end
        
        q0 = UnitQuaternion(T0);
        p0 = transl(T0);
        
        for i=1:length(s)
            qr = q0.interp(s(i));
            pr = s(i)*p0;
            T(:,:,i) = rt2tr(qr.R, pr);
        end

    else
        error('SMTB:trinterp:badarg', 'must be 2 or 3 arguments');
    end    
