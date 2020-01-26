%TREXP Matrix exponential for so(3) and se(3)
%
% For so(3)::
%
% R = TREXP(OMEGA) is the matrix exponential (3x3) of the so(3) element OMEGA that
% yields a rotation matrix (3x3). 
%
% R = TREXP(OMEGA, THETA) as above, but so(3) motion of THETA*OMEGA.
%
% R = TREXP(S, THETA) as above, but rotation of THETA about the unit vector S.
%
% R = TREXP(W) as above, but the so(3) value is expressed as a vector W
% (1x3) where W = S * THETA. Rotation by ||W|| about the vector W.
%
% For se(3)::
%
% T = TREXP(SIGMA) is the matrix exponential (4x4) of the se(3) element SIGMA that
% yields a homogeneous transformation  matrix (4x4). 
%
% T = TREXP(SIGMA, THETA) as above, but se(3) motion of SIGMA*THETA, the
% rotation part of SIGMA (4x4) must be unit norm.
%
% T = TREXP(TW) as above, but the se(3) value is expressed as a twist vector TW
% (1x6). 
%
% T = TREXP(TW, THETA) as above, but se(3) motion of TW*THETA, the
% rotation part of TW (1x6) must be unit norm.
%
% Notes::
% - Efficient closed-form solution of the matrix exponential for arguments that are
%   so(3) or se(3).
% - If THETA is given then the first argument must be a unit vector or a
%   skew-symmetric matrix from a unit vector.
% - Angle vector argument order is different to ANGVEC2R.
%
% References::
% - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p42-43.
% - Mechanics, planning and control, Park & Lynch, Cambridge, 2017.
%
% See also ANGVEC2R, TRLOG, TREXP2, SKEW, SKEWA, Twist.

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
function T = trexp(S, theta)

    if ishomog(S) || isvec(S,6)
        % input is se(3)
        
        if nargin == 1
            % twist vector 1x6 or augmented skew matrix 4x4
            if isvec(S,6)
                % it's a twist vector
                S = skewa(S);
            end
            T = expm(S);
        else
            % se(3) plus twist
            if all(size(S) == 4)
                % it's se(3) matrix
                [skw,v] = tr2rt(S);
            else
                % it's a twist vector
                v = S(1:3); v= v(:);
                skw = skew(S(4:6));
            end
            
            % use an efficient solution
            R = trexp(skw, theta);
            t = (eye(3,3)*theta + (1-cos(theta))*skw + (theta-sin(theta))*skw^2)*v;
            
            T = rt2tr(R,t);
        end         
    elseif isrot(S) || isvec(S,3)
        % input is so(3)
        
        if isrot(S)
            % input is 3x3 skew symmetric
            w = vex(S);
        elseif isvec(S)
            % input is a 3-vector
             w = S;   
        end
        
        % for a zero so(3) return unit matrix, theta not relevant
        if norm(w) < 10*eps
            T = eye(3,3);
            return;
        end
        
        if nargin == 1
            %  theta is not given, extract it
            theta = norm(w);
            w = unit(w);
        else
            % theta is given
            assert(isunit(w), 'SMTB:trexp:badarg',  'w must be a unit twist');
        end
        
        S = skew(w);
        
        T = eye(3,3) + sin(theta)*S + (1-cos(theta))*S^2;
        
    else
        error('SMTB:trexp:badarg', 'first argument must be so(3), 3-vector, se(3) or 6-vector');
    end
end
