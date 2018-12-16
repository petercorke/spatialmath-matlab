%ANGDIFF Difference of two angles
%
% ANGDIFF(TH1, TH2) is the difference between angles TH1 and TH2
% on the circle.  The result is in the interval [-pi pi).  Either or both
% arguments can be a vector:
% - If TH1 is a vector, and TH2 a scalar then return a vector where TH2 is modulo 
%   subtracted from the corresponding elements of TH1.
% - If TH1 is a scalar, and TH2 a vector then return a vector where the
%   corresponding elements of TH2 are modulo subtracted from TH1.
% - If TH1 and TH2 are vectors then return a vector whose elements are the modulo 
%   difference of the corresponding elements of TH1 and TH2.
%
% ANGDIFF(TH) as above but TH=[TH1 TH2].
%
% ANGDIFF(TH) is the equivalent angle to TH in the interval [-pi pi).
%
% Notes::
% - If TH1 and TH2 are both vectors they should have the same
%   orientation, which the output will assume.
%


% Copyright (C) 1993-2017, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for MATLAB (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.
%
% http://www.petercorke.com

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
                assert(all(size(th1) == size(th2)), 'RTB:angdiff:badarg', 'vectors must be same shape');
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
