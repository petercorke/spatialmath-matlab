%ANGDIFF Difference of two angles
%
% ANGDIFF(TH1, TH2) is the difference between angles TH1 and TH2, ie. TH1-TH2
% on the circle.  The result is in the interval [-pi pi).  Either or both
% arguments can be a vector:
% - If TH1 is a vector, and TH2 a scalar then return a vector where TH2 is modulo 
%   subtracted from the corresponding elements of TH1.
% - If TH1 is a scalar, and TH2 a vector then return a vector where the
%   corresponding elements of TH2 are modulo subtracted from TH1.
% - If TH1 and TH2 are vectors then return a vector whose elements are the modulo 
%   difference of the corresponding elements of TH1 and TH2, which must be the 
%   same length.
%
% ANGDIFF(TH) as above but TH=[TH1 TH2].
%
% ANGDIFF(TH) is the equivalent angle to the scalar TH in the interval [-pi pi).
%
% Notes::
% - The MathWorks Robotics Systems Toolbox defines a function with the same name
%   which computes TH2-TH1 rather than TH1-TH2.
% - If TH1 and TH2 are both vectors they should have the same
%   orientation, which the output will assume.
%

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

function d = angdiff(th1, th2)
  
    switch nargin
        case 1
            if length(th1) == 2
                d = th1(1) - th1(2);
            else
                d = th1;
            end
        case 2
            if length(th1) > 1 && length(th2) > 1
                % if both arguments are vectors, they must be the same
                assert(all(size(th1) == size(th2)), 'SMTB:angdiff:badarg', 'vectors must be same shape');
            end
            % th1 or th2 could be scalar
            d = th1 - th2;
    end
    
    % wrap the result into the interval [-pi pi)
    d = mod(d+pi, 2*pi) - pi;
end

% Simplistic version of the code, easy to see what it does, but slow...
%
% for very negative angles keep adding 2pi
%     while true
%         k = find(d < -pi);
%         if isempty(k)
%             break;
%         end
%         d(k) = d(k) + 2*pi;
%     end
% 
%     % for very positive angles keep subtracting 2pi
%     while true
%         k = find(d > pi);
%         if isempty(k)
%             break;
%         end
%         d(k) = d(k) - 2*pi;
%     end
