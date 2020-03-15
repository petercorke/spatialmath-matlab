%TRLOG Logarithm of SO(3) or SE(3) matrix
%
% S = trlog(R) is the matrix logarithm (3x3) of R (3x3) which is a skew
% symmetric matrix corresponding to the vector theta*w where theta is the
% rotation angle and w (3x1) is a unit-vector indicating the rotation axis.
%
% [theta,w] = trlog(R) as above but returns directly theta the rotation
% angle and w (3x1) the unit-vector indicating the rotation axis.
%
% S = trlog(T) is the matrix logarithm (4x4) of T (4x4) which has a 
% skew-symmetric upper-left 3x3 submatrix corresponding to the vector
% theta*w where theta is the rotation angle and w (3x1) is a unit-vector
% indicating the rotation axis, and a translation component.
%
% [theta,twist] = trlog(T) as above but returns directly theta the rotation
% angle and a twist vector (6x1) comprising [v w].
%
% Notes::
% - Efficient closed-form solution of the matrix logarithm for arguments that are
%   SO(3) or SE(3).
% - Special cases of rotation by odd multiples of pi are handled.
% - Angle is always in the interval [0,pi].
% - There is no Toolbox function for SO(2) or SE(2), use LOGM instead.
%
% References::
% - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p43.
% - Mechanics, planning and control, Park & Lynch, Cambridge, 2016.
%
% See also trexp, trexp2, Twist, LOGM.

%## 3d homogeneous differential

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

function [o1,o2] = trlog(T)
    
    if isa(T, 'symfun')
        T = formula(T);
    end
    
    if isrot(T)
        % deal with rotation matrix
        
        % closed form solution for matrix logarithm of a homogeneous transformation (Park & Lynch)
        % that handles the special cases
        
        % for now assumes T0 is the world frame
        
        R = T;
        
        if ~isa(R, 'sym') && abs(trace(R) - 3) < 100*eps
            % matrix is identity
            
            w = [0 0 0]';
            theta = 0;
            
        elseif ~isa(R, 'sym') && abs(trace(R) + 1) < 100*eps
            % tr R = -1
            % rotation by +/- pi, +/- 3pi etc
            
            [mx,k] = max(diag(R));
            I = eye(3,3);
            col = R(:,k) + I(:,k);
            w = col / sqrt(2*(1+mx));
            
            theta = pi;
            
%             skw = logm(R);
%             w = vex( skw );
            
        else
            % general case
            theta = acos( (trace(R)-1)/2 );
            
            skw = (R-R')/2/sin(theta);
            w = vex(skw);   % is a unit vector
            
        end
        
        if nargout <= 1
            o1 = skew(w*theta);
        elseif nargout ==2
            o1 = theta;
            o2 = w;
        end

    elseif ishomog(T)
        % SE(3) matrix
        
        [R,t] = tr2rt(T);
        
        if ~isa(R, 'sym') && abs(trace(R) - 3) < 100*eps
            % rotation matrix is identity, theta=0
            w = [0 0 0]';
            v = t;
            theta = 1;
            skw = zeros(3,3);
            
        else
            [theta, w] = trlog(R);
            skw = skew(w);
            
            Ginv = eye(3,3)/theta - skw/2 + ( 1/theta - cot(theta/2)/2 )*skw^2;
            v = Ginv * t;
        end
        
        if nargout <= 1
            o1 = [skw v; 0 0 0 0]*theta;
        elseif nargout ==2
            o1 = theta;
            o2 = [v; w];
        end
    else
        error('SMTB:trlog:badarg', 'expect SO(3) or SE(3) matrix');
    end
        
end

%     [th,w] = tr2angvec(R);
%     w = w'
%
%     d = dot(unit(w), transl(T))
%     h = d / th
%
%     q = (transl(T) - h*th*w ) * inv(eye(3,3) - R)
%
%     v =
%     rho = (eye(3,3) - R')*t / 2 / (1-cos(th));
%
%     v = cross(rho, w);
%
%     tw = [skew(unit(w)) v'; 0 0 0  0];
