%PLOTVOL Set the bounds for a 2D or 3D plot
%
% 2D plots::
%
% PLOTVOL([WX WY]) creates a new axis, and sets the bounds for a 2D plot,
% with X spanning [-WX, WX] and Y spanning [-WY,WY].
%
% PLOTVOL([XMIN XMAX YMIN YMAX]) as above but the X and Y axis limits are explicitly provided.
%
% 3D plots::
%
% PLOTVOL(W) creates a new axis, and sets the bounds for a 3D plot with X, Y and Z spanning
% the interval -W to W.  
%
% PLOTVOL([WX WY WZ]) as above with X spanning [-WX, WX], Y spanning [-WY, WY] and Z 
% spanning [-WZ, WZ].
%
% Notes::
% - The axes are labelled, grid is enabled, aspect ratio set to 1:1.
% - Hold is enabled for subsequent plots.
%
% See also: axis.

%
% PLOTVOL([XMIN XMAX YMIN YMAX ZMIN ZMAX]) as above but the X, Y and Z axis limits are
% explicitly provided.
%
% See also axis, xaxis, yaxis.

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

function plotvol(bounds)
    
    assert(any(length(bounds) == [1 2 3 4 6]), 'SMTB:plotvol:badarg', 'expecting vector of length 1, 2, 3, 4 or 6');
    
    
    clf
    axis equal
    
    switch length(bounds)
        case 1
            bounds = bounds * [-1 1 -1 1 -1 1];
        case 3
            bounds = bounds(:) * [-1 1];
            bounds = bounds';
            bounds = bounds(:)';
        case 6
            %         % 3D case
            %         set(gca, 'XTick', floor(bounds(1)):floor(bounds(2)))
            %         set(gca, 'YTick', floor(bounds(3)):floor(bounds(4)))
            %         set(gca, 'ZTick', floor(bounds(5)):floor(bounds(6)))
            bounds = bounds(:)';
        case 2
            bounds = bounds(:) * [-1 1];
            bounds = bounds';
            bounds = bounds(:)';
        case 4
            bounds = bounds(:)';
    end
    
    axis(bounds)
    if length(bounds) == 6
        view(3)
    end
    
    % common code
    grid on;
    hold on
    
    
    xyzlabel
end
