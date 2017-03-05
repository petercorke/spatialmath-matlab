%PLOT_POINT Draw a point
%
% PLOT_POINT(P, OPTIONS) adds point markers to the current plot, where P (2xN)
% and each column is the point coordinate.
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
% Examples::
%   Simple point plot
%        P = rand(2,4);
%        plot_point(P);
%
%   Plot points with markers
%        plot_point(P, '*');
%
%   Plot points with markers
%        plot_point(P, 'o', 'MarkerFaceColor', 'b');
%
%   Plot points with square markers and labels 1 to 4
%        plot_point(P, 'sequence', 's');
%
%   Plot points with circles and annotations P1 to P4
%        data = [1 2 4 8];
%        plot_point(P, 'printf', {' P%d', data}, 'o');
%
% Notes::
% - The point(s) and annotations are added to the current plot.
% - Points can be drawn in 3D axes but will always lie in the
%   xy-plane.
%
% See also PLOT, TEXT.


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

    holdon = ishold();
	hold on
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
