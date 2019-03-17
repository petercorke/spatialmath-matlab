%HOMTRANS Apply a homogeneous transformation
%
% P2 = HOMTRANS(T, P) applies the homogeneous transformation T to the points 
% stored columnwise in P.
%
% - If T is in SE(2) (3x3) and
%   - P is 2xN (2D points) they are considered Euclidean (R^2)
%   - P is 3xN (2D points) they are considered projective (P^2)
% - If T is in SE(3) (4x4) and
%   - P is 3xN (3D points) they are considered Euclidean (R^3)
%   - P is 4xN (3D points) they are considered projective (P^3)
%
% P2 and P have the same number of rows, ie. if Euclidean points are given
% then Euclidean points are returned, if projective points are given then
% projective points are returned.
%
% TP = HOMTRANS(T, T1) applies homogeneous transformation T to the 
% homogeneous transformation T1, that is TP=T*T1.  If T1 is a 3-dimensional 
% transformation then T is applied to each plane as defined by the first two 
% dimensions, ie. if T is NxN and T1 is NxNxM then the result is NxNxM.
%
% Notes::
% - If T is a homogeneous transformation defining the pose of {B} with respect to {A},
%   then the points are defined with respect to frame {B} and are transformed to be
%   with respect to frame {A}.
%
% See also E2H, H2E, RTBPose.mtimes.

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
function pt = homtrans(T, p)

    if numrows(p) == numrows(T)
        if ndims(p) == 3
            % homtrans(T1, T2)
            pt = [];
            for i=1:size(p,3)
                pt = cat(3, pt, T*p(:,:,i));
            end
        else
            % points are in projective coordinates
            pt = T * p;
        end
    elseif (numrows(T)-numrows(p)) == 1
        % points are in Euclidean coordinates, promote to homogeneous
        pt = h2e( T * e2h(p) );
    else
        error('SMTB:homtrans:badarg', 'matrices and point data do not conform')
    end
