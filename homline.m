%HOMLINE Homogeneous line from two points
%
% L = HOMLINE(X1, Y1, X2, Y2) is a vector (3x1) which describes a line in
% homogeneous form that contains the two Euclidean points (X1,Y1) and (X2,Y2).
%
% Homogeneous points X (3x1) on the line must satisfy L'*X = 0.
%
% See also PLOT_HOMLINE.

% TODO, probably should be part of a HomLine class.

function l = homline(x1, y1, x2, y2)

    l = cross([x1 y1 1], [x2 y2 1]);

    % normalize so that the result of x*l' is the pixel distance
    % from the line
    l = l / norm(l(1:2));
