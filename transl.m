%TRANSL SE(3) translational homogeneous transform
%
% Create a translational SE(3) matrix::
%
% T = TRANSL(X, Y, Z) is an SE(3) homogeneous transform (4x4) representing
% a pure translation of X, Y and Z.
%
% T = TRANSL(P) is an SE(3) homogeneous transform (4x4) representing a
% translation of P=[X,Y,Z]. P (Mx3) represents a sequence and T
% (4x4xM) is a sequence of homogeneous transforms such that T(:,:,i)
% corresponds to the i'th row of P.
%
% Extract the translational part of an SE(3) matrix::
%
% P = TRANSL(T) is the translational part of a homogeneous transform T as a
% 3-element column vector.  T (4x4xM) is a homogeneous transform
% sequence and the rows of P (Mx3) are the translational component of the
% corresponding transform in the sequence.
%
% [X,Y,Z] = TRANSL(T) is the translational part of a homogeneous transform
% T as three components.  If T (4x4xM) is a homogeneous transform sequence
% then X,Y,Z (1xM) are the translational components of the corresponding
% transform in the sequence.
%
% Notes::
% - Somewhat unusually, this function performs a function and its inverse.  An
%   historical anomaly.
%
% See also SE3.t, SE3.transl.

%## 3d homogeneous translation

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

function [t1,t2,t3] = transl(x, y, z)
    if nargin == 1
        if ishomog(x)
            if ndims(x) == 3
                % transl(T)  -> P, trajectory case
                if nargout == 1 
                    t1 = squeeze(x(1:3,4,:))';
                elseif nargout == 3
                    t1 = squeeze(x(1,4,:))';
                    t2 = squeeze(x(2,4,:))';
                    t3 = squeeze(x(3,4,:))';
                end
            else
                % transl(T)  -> P
                if nargout == 1 || nargout == 0
                    t1 = x(1:3,4);
                elseif nargout == 3
                    t1 = x(1,4);
                    t2 = x(2,4);
                    t3 = x(3,4);
                end
                    
            end
        elseif numel(x) == 3
            % transl(P) -> T
            t = x(:);
            t1 =    [eye(3)          t(:);
                0   0   0   1];
        else
            % transl(P) -> T, trajectory case
            n = size(x,1);
            t1 = repmat(eye(4,4), [1 1 n]);
            t1(1:3,4,:) = x';
        end    
    elseif nargin == 3
        % transl(x,y,z) -> T
        t = [x; y; z];
        t1 =    [ eye(3) t; 0 0 0 1];
    else
        error('SMTB:transl:badargs', 'should be 1 or 3 arguments');
    end
end

% one less function to upload for Cody/LTI assessments
function h = ishomog(tr, rtest)
    d = size(tr);
    if ndims(tr) >= 2
        h =  all(d(1:2) == [4 4]);

        if h && nargin > 1
            h = abs(det(tr(1:3,1:3)) - 1) < eps;
        end

    else
        h = false;
    end
end
