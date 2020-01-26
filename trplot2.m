%TRPLOT2 Plot a 2D coordinate frame
%
% TRPLOT2(T, OPTIONS) draws a 2D coordinate frame represented by the SE(2)
% homogeneous transform T (3x3).
%
% H = TRPLOT2(T, OPTIONS) as above but returns a handle.
%
% TRPLOT(R, OPTIONS) as above but the coordinate frame is rotated about the
% origin according to the orthonormal rotation matrix R (2x2).
%
% H = TRPLOT(R, OPTIONS) as above but returns a handle.
%
% H = TRPLOT2() creates a default frame EYE(2,2) at the origin and returns a
% handle.
%
% Animation::
%
% Firstly, create a plot and keep the the handle as per above.
%
% TRPLOT2(H, T) moves the coordinate frame described by the handle H to
% the SE(2) pose T (3x3).
%
% Options::
% 'handle',h             Update the specified handle
% 'axhandle',A           Draw in the MATLAB axes specified by the axis handle A
%
% 'color', c             The color to draw the axes, MATLAB ColorSpec
% 'axes'                 Show the MATLAB axes, box and ticks (default true)
% 'axis',A               Set dimensions of the MATLAB axes to A=[xmin xmax ymin ymax]
% 'frame',F              The frame is named {F} and the subscript on the axis labels is F.
% 'framelabel',F         The coordinate frame is named {F}, axes have no subscripts.
% 'framelabeloffset',O   Offset O=[DX DY] frame labels in units of text box height
% 'text_opts', opt       A cell array of Matlab text properties
% 'length',s             Length of the coordinate frame arms (default 1)
% 'thick',t              Thickness of lines (default 0.5)
% 'text'                 Enable display of X,Y,Z labels on the frame (default true)
% 'labels',L             Label the X,Y,Z axes with the 1st and 2nd character of the string L
% 'arrow'                Use arrows rather than line segments for the axes
% 'width', w             Width of arrow tips
% 'lefty'                Draw left-handed frame (dangerous)
%
% Examples::
%
%       trplot2(T, 'frame', 'A')
%       trplot2(T, 'frame', 'A', 'color', 'b')
%       trplot2(T1, 'frame', 'A', 'text_opts', {'FontSize', 10, 'FontWeight', 'bold'})
%
% Notes::
% - Multiple frames can be added using the HOLD command
% - When animating a coordinate frame it is best to set the axis bounds initially.
% - The 'arrow' option requires https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3
%
% See also TRPLOT.

%## 2d homogeneous graphics

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

%   'frame', name   name of the frame, used for axis subscripts and origin
%   'color', color  Matlab color specificication for the frame and annotations
%   'noaxes'        show the frame but no Matlab axes
%   'arrow'         use the contributed arrow3 function to draw the frame axes
%   'width', width  width of lines to draw if using arrow3

function hout = trplot2(T, varargin)

    if nargin == 0
        T = eye(2,2);
    end
    
    opt.color = 'b';
    opt.textcolor = [];
    opt.axes = true;
    opt.axis = [];
    opt.frame = [];
    opt.framelabel = [];
    opt.text_opts = [];
    opt.width = 1;
    opt.arrow = false;
    opt.handle = [];
    opt.axhandle = [];
    opt.length = 1;
    opt.lefty = false;
    opt.framelabeloffset = 0.2*[1 1];
    opt.handle = [];
    opt.thick = 0.5;
    opt.labels = 'XY';
    opt.text = true;


    [opt,args] = tb_optparse(opt, varargin);

    if opt.arrow && ~exist('arrow3')
        opt.arrow = false;
        warning('SMTB:trplot:badarg', 'arrow option requires arrow3 from FileExchange');
    end
    
    if isscalar(T) && ishandle(T)
        warning('SMTB:trplot2:deprecated', 'Use ''handle'' option');
        % trplot(H, T)
        opt.handle = T; T = args{1};
    end
    
    % ensure it's SE(2)
    if all(size(T) == [2 2])
        T = [T [0; 0]; 0 0 1];
    end
    
    if ~isempty(opt.handle)
        set(opt.handle, 'Matrix', se2t3(T));
        if nargout > 0
            hout = opt.handle;
        end
        return;
    end
    
    if isempty(opt.textcolor)
        opt.textcolor = opt.color;
    end
   
    if isempty(opt.text_opts)
        opt.text_opts = {};
    end
    
    % figure the dimensions of the axes, if not given
    if isempty(opt.axis)
        if all(size(T) == [3 3]) || norm(transl(T)) < eps
            c = transl2(T);
            d = 1.2;
            opt.axis = [c(1)-d c(1)+d c(2)-d c(2)+d];
        end
    end

    if ~isempty(opt.axhandle)
        hax = opt.axhandle;
        hold(hax);
        ih = ishold;
    else
        ih = ishold;
        if ~ih
            % if hold is not on, then clear the axes and set scaling
            cla
            if ~isempty(opt.axis)
                axis(opt.axis);
            end
            %axis equal
            daspect([1 1 1])
            
            if opt.axes
                xlabel( 'X');
                ylabel( 'Y');
            end
        end
        hax = gca;
        hold on
    end
    
    % create unit vectors
    o =  opt.length*[0 0 1]'; o = o(1:2);
    x1 = opt.length*[1 0 1]'; x1 = x1(1:2);
    
    if opt.lefty
        y1 = opt.length*[0 -1 1]'; y1 = y1(1:2);
    else
        y1 = opt.length*[0 1 1]'; y1 = y1(1:2);
    end
    
    opt.text_opts = [opt.text_opts, 'Color', opt.color];
    
    
    % draw the axes
    
    mstart = [o o]';
    mend = [x1 y1]';
    
    hg = hgtransform('Parent', hax);
    if opt.arrow
        % draw the 2 arrows
        S = [opt.color num2str(opt.width)];
        daspect([1 1 1]);
        ha = arrow3(mstart, mend, S);
        for h=ha'
            set(h, 'Parent', hg);
        end
    else
        for i=1:2
            plot2([mstart(i,1:2); mend(i,1:2)], 'Color', opt.color, ...
                'LineWidth', opt.thick, ...
                'Parent', hg);
        end
    end

    % label the axes
    if isempty(opt.frame)
        fmt = '%c';
    else
        fmt = sprintf('%%c_{%s}', opt.frame);
    end
    
    if opt.text
        % add the labels to each axis
        h = text(x1(1), x1(2), sprintf(fmt, opt.labels(1)), 'Parent', hg);
        if ~isempty(opt.text_opts)
            set(h, opt.text_opts{:});
        end
    end

    if opt.arrow
        set(h, 'Parent', hg);
    end

    h = text(y1(1), y1(2), sprintf(fmt, opt.labels(2)), 'Parent', hg);
    if ~isempty(opt.text_opts)
        set(h, opt.text_opts{:});
    end
    if opt.arrow
        set(h, 'Parent', hg);
    end

    if ~isempty(opt.framelabel)
        opt.frame = opt.framelabel;
    end
    % label the frame
    if ~isempty(opt.frame)
        h = text(o(1), o(2), ...
            ['\{' opt.frame '\}'], 'Parent', hg);
        set(h, 'VerticalAlignment', 'middle', ...
            'Color', opt.textcolor, ...
            'HorizontalAlignment', 'center', opt.text_opts{:});
        e = get(h, 'Extent');
        d = e(4); % use height of text box as a scale factor
        e(1:2) = e(1:2) - opt.framelabeloffset * d;
        set(h, 'Position', e(1:2));
    end
    
    if ~opt.axes
        set(gca, 'visible', 'off');
    end
    if ~ih
        hold off
    end
    
    % now place the frame in the desired pose
    set(hg, 'Matrix', se2t3(T));

    if nargout > 0
        hout = hg;
    end
end

function T3 = se2t3(T2)
    T3 = [T2(1:2,1:2) zeros(2,1) T2(1:2,3); 0 0 1 0; 0 0 0 1];
end
