%TRANSL2 SE(2) translational homogeneous transform
%
% Create a translational SE(2) matrix::
%
% T = TRANSL2(X, Y) is an SE(2) homogeneous transform (3x3) representing a
% pure translation.
%
% T = TRANSL2(P) is a homogeneous transform representing a translation or
% point P=[X,Y]. P (Mx2) represents a sequence and T (3x3xM) is a
% sequence of homogenous transforms such that T(:,:,i) corresponds to the
% i'th row of P.
%
% Extract the translational part of an SE(2) matrix::
%
% P = TRANSL2(T) is the translational part of a homogeneous transform as a
% 2-element column vector.  T (3x3xM) is a homogeneous transform
% sequence and the rows of P (Mx2) are the translational component of the
% corresponding transform in the sequence.
%
% Notes::
% - Somewhat unusually, this function performs a function and its inverse.  An
%   historical anomaly.
%
% See also SE2.t, ROT2, ISHOMOG2, TRPLOT2, TRANSL.

%## 2d homogeneous translation

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

function T = transl2(x, y)
    if nargin == 1
        if ishomog2(x)
            if ndims(x) == 3
                % transl(T)  -> P, trajectory case
                T = squeeze(x(1:2,3,:))';
            else
                % transl(T)  -> P
                T = x(1:2,3);
            end
        elseif all(size(x) == [3 3])
            T = x(1:2,3);
        elseif length(x) == 2
            % transl(P) -> T
            t = x(:);
            T =    [eye(2)          t(:);
                0   0   1];
        else
            % transl(P) -> T, trajectory case
            n = size(x,1);
            T = repmat(eye(3,3), [1 1 n]);
            T(1:2,3,:) = x';
        end    
    elseif nargin == 2
        % transl(x,y) -> T
        t = [x; y];
        T =    [ eye(2) t; 0 0 1];        
    end
end

% one less function to upload for Cody/LTI assessments

function h = ishomog2(tr, rtest)
    d = size(tr);
    if ndims(tr) >= 2
        h =  all(d(1:2) == [3 3]);

        if h && nargin > 1
            h = abs(det(tr(1:2,1:2)) - 1) < eps;
        end
    else
        h = false;
    end
end
