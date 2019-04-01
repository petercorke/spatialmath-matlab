%RPY2TR Roll-pitch-yaw angles to SE(3) homogeneous transform
%
% T = RPY2TR(ROLL, PITCH, YAW, OPTIONS) is an SE(3) homogeneous
% transformation matrix (4x4) with zero translation and rotation equivalent
% to the specified roll, pitch, yaw angles angles. These correspond to
% rotations about the Z, Y, X axes respectively. If ROLL, PITCH, YAW are
% column vectors (Nx1) then they are assumed to represent a trajectory and
% R is a three-dimensional matrix (4x4xN), where the last index corresponds
% to rows of ROLL, PITCH, YAW.
%
% T = RPY2TR(RPY, OPTIONS) as above but the roll, pitch, yaw angles are
% taken from the vector (1x3) RPY=[ROLL,PITCH,YAW]. If RPY is a matrix
% (Nx3) then R is a three-dimensional matrix (4x4xN), where the last index
% corresponds to rows of RPY which are assumed to be
% ROLL,PITCH,YAW].
%
% Options::
%  'deg'      Compute angles in degrees (radians default)
%  'xyz'      Rotations about X, Y, Z axes (for a robot gripper)
%  'zyx'      Rotations about Z, Y, X axes (for a mobile robot, default)
%  'yxz'      Rotations about Y, X, Z axes (for a camera)
%  'arm'      Rotations about X, Y, Z axes (for a robot arm)
%  'vehicle'  Rotations about Z, Y, X axes (for a mobile robot)
%  'camera'   Rotations about Y, X, Z axes (for a camera)
%
% Note::
% - Toolbox rel 8-9 has the reverse angle sequence as default.
% - ZYX order is appropriate for vehicles with direction of travel in the X
%   direction.  XYZ order is appropriate if direction of travel is in the Z
%   direction.
% - 'arm', 'vehicle', 'camera' are synonyms for 'xyz', 'zyx' and 'yxz'
%   respectively.
%
% See also TR2RPY, RPY2R, EUL2TR.

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

function T = rpy2tr(roll, varargin)

    R = rpy2r(roll, varargin{:});
    T = r2t(R);
