%TR2ANGVEC Convert rotation matrix to angle-vector form
%
% [THETA,V] = TR2ANGVEC(R, OPTIONS) is rotation expressed in terms of an
% angle THETA (1x1) about the axis V (1x3) equivalent to the orthonormal rotation
% matrix R (3x3).
%
% [THETA,V] = TR2ANGVEC(T, OPTIONS) as above but uses the rotational part of the
% homogeneous transform T (4x4).
%
% If R (3x3xK) or T (4x4xK) represent a sequence then THETA (Kx1)is a vector 
% of angles for corresponding elements of the sequence and V (Kx3) are the 
% corresponding axes, one per row.
%
% Options::
% 'deg'   Return angle in degrees (default radians)
%
% Notes::
% - For an identity rotation matrix both THETA and V are set to zero.
% - The rotation angle is always in the interval [0 pi], negative rotation
%   is handled by inverting the direction of the rotation axis.
% - If no output arguments are specified the result is displayed.
%
% See also ANGVEC2R, ANGVEC2TR, TRLOG.

%## 3d homogeneous rotation

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

function [theta_, n_] = tr2angvec(R, varargin)

    assert(ishomog(R) || isrot(R), 'SMTB:tr2angvec:badarg', 'argument must be SO(3) or SE(3)');

    opt.deg = false;
    opt = tb_optparse(opt, varargin);
    
    % get the rotation submatrix(s)
    if ~isrot(R)
        R = t2r(R);
    end
    
    if size(R,3) > 1
        theta = zeros(size(R,3),1);
        v = zeros(size(R,3),3);
    end
    
    for i=1:size(R,3)  % for each rotation matrix in the sequence
        
        % There are a few ways to do this:
        %
        % 1.
        %
        % e = 0.5*vex(R - R');  % R-R' is skew symmetric
        % theta = asin(norm(e));
        % n = unit(e);
        %
        %  but this fails for rotations > pi/2
        %
        % 2.
        %
        % e = vex(logm(R));
        % theta = norm(e);
        % n = unit(e);
        %
        %  elegant, but 40x slower than using eig
        %
        % 3.
        %
        % Use eigenvectors, get angle from trace which is defined over -pi to
        % pi.  Don't use eigenvalues since they only give angles -pi/2 to pi/2.
        %
        % 4.
        %
        % Take the log of the rotation matrix
        
        Ri = R(:,:,i);
        
        % check the determinant
        assert( abs(det(Ri)-1) < 10*eps, 'SMTB:tr2angvec:badarg', 'matrix is not orthonormal');
        
        [th,v] = trlog(Ri);
        theta(i) = th;
        n(i,:) = v;
        
        if opt.deg
            theta(i) = theta(i) * 180/pi;
            units = 'deg';
        else
            units = 'rad';
        end
        
        if nargout == 0
            % if no output arguments display the angle and vector
            fprintf('Rotation: %f %s x [%f %f %f]\n', theta(i), units, n(i,:));
        end
    end
    
    if nargout == 1
        theta_ = theta;
    elseif nargout == 2
        theta_ = theta;
        n_ = n;
    end
