%TRSCALE Homogeneous transformation for pure scale
%
% T = TRSCALE(S) is a homogeneous transform (4x4) corresponding to a pure
% scale change.  If S is a scalar the same scale factor is used for x,y,z,
% else it can be a 3-vector specifying scale in the x-, y- and
% z-directions.
%
% Note::
% - This matrix does not belong to SE(3) and if compounded with
%   any SE(3) matrix the result will not be in SE(3).
%

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

function t = trscale(sx, sy, sz)

    if length(sx) > 1
        s = sx;
    else
        if nargin == 1
            s = [sx sx sx];
        else
            s = [sx sy sz];
        end
    end
    t = r2t(diag(s));
