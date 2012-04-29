%ANGDIFF Difference of two angles
%
% D = ANGDIFF(TH1, TH2) returns the difference between angles TH1 and TH2 on
% the circle.  The result is in the interval [-pi pi).  If TH1 is a column 
% vector, and TH2 a scalar then return a column vector where TH2 is modulo 
% subtracted from the corresponding elements of TH1.
%
% D = ANGDIFF(TH) returns the equivalent angle to TH in the interval [-pi pi).
%
% Return the equivalent angle in the interval [-pi pi).

function d = angdiff(th1, th2)

    if nargin < 2
% THIS IS A BAD IDEA, WHERE IS IT USED?
%         if length(th1) > 1
%             d = th1(1) - th1(2);
%         else
%             d = th1;
%         end
        d = th1;
    else
        d = th1 - th2;
    end

    
    d = mod(d+pi, 2*pi) - pi;

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
