%T2R Rotational submatrix
%
% R = T2R(T) is the orthonormal rotation matrix component of homogeneous
% transformation matrix T.  Works for T in SE(2) or SE(3)
% - If T is 4x4, then R is 3x3.
% - If T is 3x3, then R is 2x2.
%
% Notes::
% - For a homogeneous transform sequence (KxKxN) returns a rotation matrix
%   sequence (K-1xK-1xN).
% - The validity of rotational part is not checked
%
% See also R2T, TR2RT, RT2TR.

%## 2d 3d homogeneous

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

function R = t2r(T)
    
    assert(isfloat(T), 'SMTB:t2r:badarg', 'expecting real matrix argument');
    
    % check dimensions: T is SE(2) or SE(3)
    d = size(T);
    assert(d(1) == d(2), 'SMTB:t2r:badarg', 'matrix must be square');
    assert(any(d(1) == [3 4]), 'SMTB:t2r:badarg', 'argument is not a homogeneous transform (sequence)');
    
    n = d(1);     % works for SE(2) or SE(3)
    
    if numel(d) == 2
        % single matrix case
        R = T(1:n-1,1:n-1);
    else
        %  matrix sequence case
        R = zeros(n-1,n-1,d(3));
        for i=1:d(3)
            R(:,:,i) = T(1:n-1,1:n-1,i);
        end
    end
end
