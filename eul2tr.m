%EUL2TR Convert Euler angles to homogeneous transform
%
% T = EUL2TR(PHI, THETA, PSI, OPTIONS) is an SE(3) homogeneous
% transformation matrix (4x4) with zero translation and rotation equivalent
% to the specified Euler angles. These correspond to rotations about the Z,
% Y, Z axes respectively. If PHI, THETA, PSI are column vectors (Nx1) then
% they are assumed to represent a trajectory and R is a three-dimensional
% matrix (4x4xN), where the last index corresponds to rows of PHI, THETA,
% PSI.
%
% R = EUL2R(EUL, OPTIONS) as above but the Euler angles are taken from the
% vector (1x3)  EUL = [PHI THETA PSI]. If EUL is a matrix (Nx3) then R is a
% three-dimensional matrix (4x4xN), where the last index corresponds to
% rows of RPY which are assumed to be [PHI,THETA,PSI].
%
% Options::
%  'deg'      Angles given in degrees (radians default)
%
% Note::
% - The vectors PHI, THETA, PSI must be of the same length.
% - The translational part is zero.
%
% See also EUL2R, RPY2TR, TR2EUL, SE3.eul.

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

function T = eul2tr(phi, varargin)

    R = eul2r(phi, varargin{:});
    T = r2t(R);
