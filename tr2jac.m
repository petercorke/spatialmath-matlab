%TR2JAC Jacobian for differential motion
%
% J = TR2JAC(TAB) is a Jacobian matrix (6x6) that maps spatial velocity or
% differential motion from frame {A} to frame {B} where the pose of {B}
% relative to {A} is represented by the homogeneous transform TAB (4x4).
%
% J = TR2JAC(TAB, 'samebody') is a Jacobian matrix (6x6) that maps spatial
% velocity or differential motion from frame {A} to frame {B} where both
% are attached to the same moving body.  The pose of {B} relative to {A} is
% represented by the homogeneous transform TAB (4x4).
%
% See also WTRANS, TR2DELTA, DELTA2TR, SE3.velxform.

%## 3d homogeneous differential

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

function J = tr2jac(T, varargin)
    
    opt.samebody = false;
    
    opt = tb_optparse(opt, varargin);
		
    R = t2r(T);
    
    if opt.samebody
        J = [
            R'              (skew(transl(T))*R)'
            zeros(3,3)      R'
            ];
    else
        J = [
            R'              zeros(3,3)
            zeros(3,3)      R'
            ];
    end
