%VEXA Convert augmented skew-symmetric matrix to vector
%
% V = VEXA(S) is the vector which has the corresponding augmented skew-symmetric 
% matrix S.  
%
% In the case that S (3x3) = 
%
%               |  0  -v3  v1 |
%               | v3    0  v2 |
%               |  0    0   0 |
%
% then V = [v1; v2; v3].  In the case that S (6x6) = 
%
%
%               |  0  -v6   v5  v1 |
%               | v6    0  -v4  v2 |
%               |-v5   v4    0  v3 |
%               |  0    0    0   0 |
%
% then V = [v1; v2; v3; v4; v5; v6].
%
% Notes::
% - This is the inverse of the function SKEWA().
% - The matrices are the generator matrices for se(2) and se(3).  The elements
%   comprise the equivalent twist vector.
%
% References::
% - Robotics, Vision & Control: Second Edition, Chap 2,
%   P. Corke, Springer 2016.
%
% See also SKEWA, VEX, Twist.

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

function s = vexa(Omega)

    if all(size(Omega) == [4 4])
        s = [transl(Omega); vex(Omega(1:3,1:3))];
    elseif all(size(Omega) == [3 3 ])
        s = [transl2(Omega); vex(Omega(1:2,1:2))];
    else
        error('SMTB:vexa:badarg', 'argument must be a 3x3 or 4x4 matrix');
    end
    
    
end
