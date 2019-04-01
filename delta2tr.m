%DELTA2TR Convert differential motion  to SE(3) homogeneous transform
%
% T = DELTA2TR(D) is a homogeneous transform (4x4) representing differential 
% motion D (6x1). 
%
% The vector D=(dx, dy, dz, dRx, dRy, dRz) represents infinitessimal translation
% and rotation, and is an approximation to the instantaneous spatial velocity 
% multiplied by time step.
%
% Reference::
% - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p67.
%
% See also tr2delta, SE3.delta.

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

function delta = delta2tr(d)
    d = d(:);
    delta = eye(4,4) + [skew(d(4:6)) d(1:3); 0 0 0 0];
