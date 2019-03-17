%PLOT_POINT Draw a point
%
% PLOT_POINT(P, OPTIONS) adds point markers and optional annotation text
% to the current plot, where P (2xN) and each column is a point coordinate.
%
% H = PLOT_POINT(...) as above but return handles to the points.
%
% Options::
%  'textcolor', colspec     Specify color of text
%  'textsize', size         Specify size of text
%  'bold'                   Text in bold font.
%  'printf', {fmt, data}    Label points according to printf format
%                           string and corresponding element of data
%  'sequence'               Label points sequentially
%  'label',L                Label for point
%
% Additional options to PLOT can be used: 
% - standard MATLAB LineStyle such as 'r' or 'b---'
% - any MATLAB LineProperty options can be given such as 'LineWidth', 2.
%
% Notes::
% - The point(s) and annotations are added to the current plot.
% - Points can be drawn in 3D axes but will always lie in the
%   xy-plane.
% - Handles are to the points but not the annotations.
%
% Examples::
%   Simple point plot with two markers
%        P = rand(2,4);
%        plot_point(P);
%
%   Plot points with markers
%        plot_point(P, '*');
%
%   Plot points with solid blue circular markers
%        plot_point(P, 'bo', 'MarkerFaceColor', 'b');
%
%   Plot points with square markers and labelled 1 to 4
%        plot_point(P, 'sequence', 's');
%
%   Plot points with circles and labelled P1, P2, P4 and P8
%        data = [1 2 4 8];
%        plot_point(P, 'printf', {' P%d', data}, 'o');
%

%
% See also PLOT_SPHERE, PLOT, TEXT.

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

function ho = plot_point(p, varargin)

    opt.textcolor = 'k';
    opt.textsize = 12;
    opt.printf = [];
    opt.sequence = false;
    opt.bold = false;
    opt.label = [];
    opt.solid = false;
    [opt,arglist,ls] = tb_optparse(opt, varargin);

    % label is a cell array, one per point (column)
    if ~isempty(opt.label) && numcols(p) == 1
        % if single point, convert single label to a cell array
        opt.label = {opt.label};
    end
    
    % default marker style
    if isempty(ls)
        ls = {'bs'};    % blue square
    end

    % add stuff to pull .u and .v out of a vector of objects
    if ~isnumeric(p) && any(strcmp('u_', properties(p)))
        % p is an object with u_ and v_ properties
        p = [[p.u_]; [p.v_]];
    end

    if numrows(p) == 3
        error('p must have 2 rows, only supports 2D plotting')
    end
    holdon = ishold();
	hold on
    h = zeros(1,numcols(p));
    
	for i=1:numcols(p)
        if opt.solid 
            arglist = [ 'MarkerFaceColor', ls{1}(1), arglist];
        end
		h(i) = plot(p(1,i), p(2,i), ls{:}, arglist{:});
        if opt.sequence
            show(p(:,i), '%d', i, opt);
        end

        if ~isempty(opt.label)
            show(p(:,i), opt.label{i}, [], opt);
        elseif ~isempty(opt.printf)
            show(p(:,i), opt.printf{1}, opt.printf{2}(i), opt);
        end

	end
    if ~holdon
        hold off
    end
    figure(gcf)
    if nargout > 0
        ho = h;
    end
end

function show(p, fmt, val, opt)
    if opt.bold
        fw = 'bold';
    else
        fw = 'normal';
    end
    text(p(1), p(2), sprintf([' ' fmt], val), ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontUnits', 'pixels', ...
        'FontSize', opt.textsize, ...
        'FontWeight', fw, ...
        'Color', opt.textcolor);
end
