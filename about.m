%ABOUT Compact display of variable type
%
% ABOUT(X) displays a compact line that describes the class and dimensions of
% X.
%
% ABOUT X  as above but this is the command rather than functional form.
%
% Examples::
%         >> a=1;
%         >> about a
%         a [double] : 1x1 (8 bytes)
%
%         >> a = rand(5,7);
%         >> about a
%         a [double] : 5x7 (280 bytes)
%
% See also WHOS.

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
function about(varargin)
    for i=1:nargin
        var = varargin{i};
    if ischar(var)
        % invoked without parentheses
        w = evalin('caller', sprintf('whos(''%s'')', var));
        varname = var;
    else
        w = whos('var');
        varname = inputname(1);
    end
    
    if isempty(w)
        error('SMTB:about', ['cant find variable ' var])
    end
    ss = sprintf('%d', w.size(1));
    for i=2:length(w.size)
        ss = strcat(ss, sprintf('x%d', w.size(i)));
    end
    
    % build a string to show if complex or not
    if w.complex
        cmplx = '+complex';
    else
        cmplx = '';
    end
    
    % build a string to show size in convenient format
    suffix = {'bytes', 'kB', 'MB', 'GB', 'TB'};
    sz = w.bytes;
    for i=1:numel(suffix)
        if sz/1000 < 1
            break;
        end
        sz = sz/1000;
    end
    
    if i==1
        size = sprintf('%d %s', sz, suffix{i});
    else
        size = sprintf('%.1f %s', sz, suffix{i});
    end
    
    % now display the info
    fprintf('%s [%s%s] : %s (%s)\n', ...
        varname, w.class, cmplx, ss, size);

    end
