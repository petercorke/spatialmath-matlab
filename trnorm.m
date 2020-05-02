%TRNORM Normalize an SO(3) or SE(3) matrix
%
% TRNORM(R) is guaranteed to be a proper orthogonal matrix rotation
% matrix (3x3) which is "close" to the input matrix R (3x3). If R
% = [N,O,A] the O and A vectors are made unit length and the normal vector
% is formed from N = O x A, and then we ensure that O and A are orthogonal
% by O = A x N.
%
% TRNORM(T) as above but the rotational submatrix of the homogeneous
% transformation T (4x4) is normalised while the translational part is
% unchanged.
%
% If R (3x3xK) or T (4x4xK) representing a sequence then the normalisation
% is performed on each of the K planes.
%
% Notes::
% - Only the direction of A (the z-axis) is unchanged.
% - Used to prevent finite word length arithmetic causing transforms to 
%   become `unnormalized'.
% - There is no Toolbox function for SO(2) or SE(2).
%
% See also OA2TR, SO3.trnorm, SE3.trnorm.

%## 3d homogeneous 

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

function TR = trnorm(T)

    assert(ishomog(T) || isrot(T), 'SMTB:trnorm:badarg', 'expecting 3x3xN or 4x4xN hom xform');
    
    if ndims(T) == 3
        % recurse for transform sequence
        n = size(T);
        TR = zeros(n);
        for i=1:n(3)
            TR(:,:,i) = trnorm(T(:,:,i));
        end
        return
    end
    
    o = T(1:3,2); a = T(1:3,3);
    n = cross(o, a);         % N = O x A
    o = cross(a, n);         % O = A x N
    R = [unit(n) unit(o) unit(a)];
    
    if ishomog(T)
        TR = rt2tr( R, T(1:3,4) );
    elseif isrot(T)
        TR = R;
    end

