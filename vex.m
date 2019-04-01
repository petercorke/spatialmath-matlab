%VEX Convert skew-symmetric matrix to vector
%
% V = VEX(S) is the vector which has the corresponding skew-symmetric 
% matrix S.  
%
% In the case that S (2x2) =
%
%           | 0  -v |
%           | v   0 |
%
% then V = [v].  In the case that S (3x3) =
%
%           |  0  -vz   vy |
%           | vz    0  -vx |
%           |-vy   vx    0 |
%
% then V = [vx; vy; vz].
%
% Notes::
% - This is the inverse of the function SKEW().
% - Only rudimentary checking (zero diagonal) is done to ensure that the 
%   matrix is actually skew-symmetric.
% - The function takes the mean of the two elements that correspond to 
%   each unique element of the matrix.
% - The matrices are the generator matrices for so(2) and so(3).
%
% References::
% - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p25+43.
%
% See also SKEW, VEXA.

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

function v = vex(S)
%     if trace(abs(S)) > 10*eps
%         error('SMTB:vex:badarg', 'argument is not skew symmetric tr=%g', trace(abs(S)));
%     end
    if all(size(S) == [3 3])
        v = 0.5*[S(3,2)-S(2,3); S(1,3)-S(3,1); S(2,1)-S(1,2)];
    elseif all(size(S) == [2 2])
        v = 0.5*(S(2,1)-S(1,2));
    else
        error('SMTB:vex:badarg', 'argument must be a 2x2 or 3x3 matrix');
    end
