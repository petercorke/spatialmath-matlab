%RT2TR Convert rotation and translation to homogeneous transform
%
% TR = RT2TR(R, t) is a homogeneous transformation matrix (N+1xN+1) formed
% from an orthonormal rotation matrix R (NxN) and a translation vector t
% (Nx1).  Works for R in SO(2) or SO(3):
%  - If R is 2x2 and t is 2x1, then TR is 3x3
%  - If R is 3x3 and t is 3x1, then TR is 4x4
%
% For a sequence R (NxNxK) and t (NxK) results in a transform sequence (N+1xN+1xK).
%
% Notes::
% - The validity of R is not checked
%
% See also T2R, R2T, TR2RT.

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

function T = rt2tr(R, t)

    assert( size(R,1) == size(R,2), 'SMTB:rt2tr:badarg', 'R must be square');
    
    n = size(R,2);
    B = [zeros(1,n) 1];
    
    if size(R,3) > 1
        % vector case
        assert( size(R,1) == size(t,1), 'SMTB:rt2tr:badarg', 'R and t must have conforming dimensions')
        assert( size(R,3) == size(t,2), 'SMTB:rt2tr:badarg', 'For sequence size(R,3) must equal size(t,2)');
        
        T = zeros(n+1,n+1,size(R,3));
        for i=1:size(R,3)
            T(:,:,i) = [R(:,:,i) t(:,i); B];
        end
    else
        % scalar case
        assert( isvec(t, size(R,1)), 'SMTB:rt2tr:badarg', 'R and t must have conforming dimensions')
        T = [R t(:); B];
    end
end

