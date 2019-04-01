%TREXP2 Matrix exponential for so(2) and se(2)
%
% SO(2)::
%
% R = TREXP2(OMEGA) is the matrix exponential (2x2) of the so(2) element OMEGA that
% yields a rotation matrix (2x2).
%
% R = TREXP2(THETA) as above, but rotation by THETA (1x1). 
%
% SE(2)::
%
% T = TREXP2(SIGMA) is the matrix exponential (3x3) of the se(2) element
% SIGMA that yields a homogeneous transformation  matrix (3x3).
%
% T = TREXP2(SIGMA, THETA) as above, but se(2) rotation of SIGMA*THETA, the
% rotation part of SIGMA (3x3) must be unit norm.
%
% T = TREXP2(TW) as above, but the se(2) value is expressed as a vector TW
% (1x3).
%
% T = TREXP(TW, THETA) as above, but se(2) rotation of TW*THETA, the
% rotation part of TW must be unit norm.
%
% Notes::
% - Efficient closed-form solution of the matrix exponential for arguments that are
%   so(2) or se(2).
% - If THETA is given then the first argument must be a unit vector or a
%   skew-symmetric matrix from a unit vector.
%
% References::
% - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p25-26.
% - Mechanics, planning and control, Park & Lynch, Cambridge, 2017.
%
% See also TREXP, SKEW, SKEWA, Twist.

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

function T = trexp2(S, theta)

    
    if ishomog2(S) || isvec(S,3)
        % input is se(2)
        if nargin == 1
            if isvec(S,3)
                S = S(:);
                S = [skew(S(3)) S(1:2); 0 0 0];
            end
            T = expm(S);
        else
            if ishomog2(S)
                v = S(1:2,3);
                skw = S(1:2,1:2);
            else
                v = S(1:2)';
                skw = skew(S(3));
            end
            
            R = trexp2(skw, theta);
            
            t = (eye(2,2)*theta + (1-cos(theta))*skw + (theta-sin(theta))*skw^2)*v;
            
            T = rt2tr(R,t);
            %T = expm([S v; 0 0 0]*theta);
        end
    elseif isrot2(S) || isvec(S,1)
        % input is so(2)
        if isrot2(S)
            % input is 2x2 skew symmetric
            w = vex(S);
        elseif isvec(S,1)
            % input is a 1-vector
            w = S;
        end
        
        if nargin == 1
            %  theta is not given, extract it
            if norm(w) < 10*eps
                T = eye(2,2);
                return;
            end
            theta = norm(w);
            S = skew(unit(w));
        else
            if theta < 10*eps
                T = eye(2,2);
                return;
            end
            if ~isunit(w)
                error('SMTB:trexp: angular velocity must be a unit vector');
            end
            S = skew(w);
        end
        
        T = eye(2,2) + sin(theta)*S + (1-cos(theta))*S^2;
        
    else
        error('SMTB:trexp2:badarg', 'first argument must be so(2), 3-vector, se(2) or 6-vector');end
end
