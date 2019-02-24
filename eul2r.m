%EUL2R Convert Euler angles to rotation matrix
%
% R = EUL2R(PHI, THETA, PSI, OPTIONS) is an SO(3) orthonornal rotation
% matrix (3x3) equivalent to the specified Euler angles.  These correspond
% to rotations about the Z, Y, Z axes respectively. If PHI, THETA, PSI are
% column vectors (Nx1) then they are assumed to represent a trajectory and
% R is a three-dimensional matrix (3x3xN), where the last index corresponds
% to rows of PHI, THETA, PSI.
%
% R = EUL2R(EUL, OPTIONS) as above but the Euler angles are taken from the
% vector (1x3)  EUL = [PHI THETA PSI]. If EUL is a matrix (Nx3) then R is a
% three-dimensional matrix (3x3xN), where the last index corresponds to
% rows of RPY which are assumed to be [PHI,THETA,PSI].
%
% Options::
%  'deg'      Angles given in degrees (radians default)
%
% Note::
% - The vectors PHI, THETA, PSI must be of the same length.
%
% See also EUL2TR, RPY2TR, TR2EUL, SO3.eul.

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

function R = eul2r(phi, varargin)
    opt.deg = false;
    [opt,args] = tb_optparse(opt, varargin);

    % unpack the arguments
    if numcols(phi) == 3
		theta = phi(:,2);
		psi = phi(:,3);
		phi = phi(:,1);
	elseif nargin >= 3
        theta = args{1};
        psi = args{2};
    else
        error('SMTB:eul2r:badarg', 'expecting 3 inputs, 3-vector or 3-column matrix')
    end

    % optionally convert from degrees
    if opt.deg
        d2r = pi/180.0;
        phi = phi * d2r;
        theta = theta * d2r;
        psi = psi * d2r;
    end

    if numrows(phi) == 1
        R = rotz(phi) * roty(theta) * rotz(psi);
    else
        R = zeros(3,3,numrows(phi));
        for i=1:numrows(phi)
            R(:,:,i) = rotz(phi(i)) * roty(theta(i)) * rotz(psi(i));
        end
    end
