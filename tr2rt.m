%TR2RT Convert homogeneous transform to rotation and translation 
%
% [R,t] = TR2RT(TR) splits a homogeneous transformation matrix (NxN) into an 
% orthonormal rotation matrix R (MxM) and a translation vector t (Mx1), where
% N=M+1.
%
% Works for TR in SE(2) or SE(3)
%  - If TR is 4x4, then R is 3x3 and T is 3x1.
%  - If TR is 3x3, then R is 2x2 and T is 2x1.
%
% A homogeneous transform sequence TR (NxNxK) is split into rotation matrix 
% sequence R (MxMxK) and a translation sequence t (KxM).
%
%
% Notes::
% - The validity of R is not checked.
%
% See also RT2TR, R2T, T2R.

%## 2d 3d homogeneous rotation translation

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

function [R,t] = tr2rt(T)
    assert(size(T,1) == size(T,2), 'T must be square');

    n = size(T,2)-1;

    if size(T,3) > 1
        R = zeros(n,n,size(T,3));
        t = zeros(size(T,3), n);
        for i=1:size(T,3)
            R(1:n,1:n,i) = T(1:n,1:n,i);
            t(i,:) = T(1:n,n+1,i)';
        end
    else
        R = T(1:n,1:n);
        t = T(1:n,n+1);
    end
