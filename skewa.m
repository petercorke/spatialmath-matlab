%SKEWA Create augmented skew-symmetric matrix
%
% S = SKEWA(V) is an augmented skew-symmetric matrix formed from V.
% 
% If V (1x3) then S =
%
%           |  0  -v3  v1 |
%           | v3    0  v2 |
%           |  0    0   0 |
%
% and if V (1x6) then S =
%
%           |  0  -v6   v5  v1 |
%           | v6    0  -v4  v2 |
%           |-v5   v4    0  v3 |
%           |  0    0    0   0 |
%
% Notes::
% - This is the inverse of the function VEXA().
% - These are the generator matrices for the Lie algebras se(2) and se(3).
% - Map twist vectors in 2D and 3D space to se(2) and se(3).
%
% References::
% - Robotics, Vision & Control: Second Edition, Chap 2,
%   P. Corke, Springer 2016.
%
% See also SKEW, VEX, Twist.

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

function Omega = skewa(s)
    s  = s(:);
    switch length(s)
        case 3
            Omega = [skew(s(3)) s(1:2); 0 0 0];
            
        case 6
            Omega = [skew(s(4:6)) s(1:3); 0 0 0 0];
            
        otherwise
            error('SMTB:skewa:badarg', 'expecting a 3- or 6-vector');
    end
end
