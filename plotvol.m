%PLOTVOL Set the bounds for a 2D or 3D plot
%
% PLOTVOL(W) creates a new axis, and sets the bounds for a 2D plot with X and Y spanning
% the interval -W to W.  The axes are labelled, grid is enabled, aspect ratio set to 1:1,
% and hold is enabled for subsequent plots.
%
% PLOTVOL([XMIN XMAX YMIN YMAX]) as above but the X and Y axis limits are explicitly provided.
%
% PLOTVOL([XMIN XMAX YMIN YMAX ZMIN ZMAX]) as above but the X, Y and Z axis limits are 
% explicitly provided.
%
% See also axis, xaxis, yaxis.

function plotvol(bounds)

    if ~any(length(bounds) == [1 4 6])
        error('RTB:plotvol:badarg', 'expecting vector of length 1, 4 or 6');
    end
    if length(bounds) == 1
        bounds = bounds * [-1 1 -1 1];
    end
    
    axis equal
    axis(bounds);

    if length(bounds) == 4
        % 2D case
        
        set(gca, 'XTick', floor(bounds(1)):floor(bounds(2)))
        set(gca, 'YTick', floor(bounds(3)):floor(bounds(4)))
        
    elseif length(bounds) == 6
        % 3D case
        set(gca, 'XTick', floor(bounds(1)):floor(bounds(2)))
        set(gca, 'YTick', floor(bounds(3)):floor(bounds(4)))
        set(gca, 'ZTick', floor(bounds(5)):floor(bounds(6)))
        
    end
    
    % common code
    grid on;
    hold on


    xyzlabel
