%ANGVEC2R Convert angle and vector orientation to a rotation matrix
%
% R = ANGVEC2R(THETA, V) is an orthonormal rotation matrix (3x3)
% equivalent to a rotation of THETA about the vector V.
%
% Notes::
% - Uses Rodrigues' formula
% - If THETA == 0 then return identity matrix and ignore V.
% - If THETA ~= 0 then V must have a finite length.
%
% See also angvec2tr, eul2r, rpy2r, tr2angvec, trexp, SO3.angvec.

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
function R = angvec2r(theta, v)

    if nargin < 2 || ~isscalar(theta) || ~isvec(v)
        error('SMTB:angvec2r:badarg', 'bad arguments');
    end
    if ~isa(v, 'sym') && norm(v) < 10*eps
        if (abs(theta) > 0)
            error('SMTB:angvec2r:badarg', 'norm of direction is zero');
        else
            R = eye(3,3);
            return;
        end
    end
   
    % Rodrigue's equation
    
    sk = skew( unit(v) );
    R = eye(3,3) + sin(theta)*sk + (1-cos(theta))*sk^2;
end
