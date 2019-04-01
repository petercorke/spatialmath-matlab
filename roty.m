%ROTY SO(3) rotation about Y axis
%
% R = ROTY(THETA) is an SO(3) rotation matrix (3x3) representing a rotation of THETA 
% radians about the y-axis.
%
% R = ROTY(THETA, 'deg') as above but THETA is in degrees.
%
% See also TROTY, ROTX, ROTZ, ANGVEC2R, ROT2, SO3.Ry.

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

function R = roty(t, deg)
    
    assert((isreal(t) & isscalar(t)) | isa(t, 'sym'), ...
        'SMTB:roty:badarg', 'theta must be a real scalar');
    
    if nargin > 1 && strcmp(deg, 'deg')
        t = t *pi/180;
    end
    
    ct = cos(t);
    st = sin(t);
    
    % make almost zero elements exactly zero
    if ~isa(t, 'sym')
        if abs(st) < eps
            st = 0;
        end
        if abs(ct) < eps
            ct = 0;
        end
    end
    
    % create the rotation matrix
    R = [
        ct  0   st
        0   1   0
       -st  0   ct
       ];
end
