%OA2TR Convert orientation and approach vectors to homogeneous transformation
%
% T = OA2TR(O, A) is an SE(3) homogeneous tranformation (4x4) for the
% specified orientation and approach vectors (3x1) formed from 3 vectors
% such that R = [N O A] and N = O x A.
%
% Notes::
% - The rotation submatrix is guaranteed to be orthonormal so long as O and A 
%   are not parallel.
% - The vectors O and A are parallel to the Y- and Z-axes of the coordinate
%   frame respectively.
% - The translational part is zero.
%
% References::
% - Robot manipulators: mathematics, programming and control
%   Richard Paul, MIT Press, 1981.
%
% See also RPY2TR, EUL2TR, OA2R, SE3.oa.

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

function r = oa2tr(o, a)
    assert( nargin >= 2 && isvec(o) && isvec(a), 'SMTB:oa2tr:badarg', 'bad arguments');
    
	n = cross(o, a);
    o = cross(a, n);
	r = [unit(n(:)) unit(o(:)) unit(a(:)) zeros(3,1); 0 0 0 1];
