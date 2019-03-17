%EUL2JAC Euler angle rate Jacobian
%
% J = EUL2JAC(PHI, THETA, PSI) is a Jacobian matrix (3x3) that maps ZYZ Euler angle rates to 
% angular velocity at the operating point specified by the Euler angles PHI, THETA, PSI.
%
% J = EUL2JAC(EUL)  as above but the Euler angles are passed as a vector EUL=[PHI, THETA, PSI]. 
%
% Notes::
% - Used in the creation of an analytical Jacobian.
% - Angles in radians, rates in radians/sec.
%
% Reference::
% - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p232-3.
%
% See also rpy2jac, eul2r, SerialLink.jacobe.

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

function J = eul2jac(phi, theta, psi)

    if length(phi) == 3
        % eul2jac([phi theta psi])
        theta = phi(2);
        psi = phi(3);
        phi = phi(1);
    elseif nargin ~= 3
        error('SMTB:eul2jac:badarg', 'bad arguments');
    end
J = [ 0, -sin(phi), cos(phi)*sin(theta)
      0,  cos(phi), sin(phi)*sin(theta)
      1,        0,           cos(theta) ];  
        
