%TRPRINT2 Compact display of SE(2) homogeneous transformation
%
% TRPRINT2(T, OPTIONS) displays the homogoneous transform (3x3) in a compact 
% single-line format.  If T is a homogeneous transform sequence then each 
% element is printed on a separate line.
%
% TRPRINT2(R, OPTIONS) as above but displays the SO(2) rotation matrix (3x3).
%
% S = TRPRINT2(T, OPTIONS) as above but returns the string.
%
% TRPRINT2 T  is the command line form of above, and displays in RPY format.
%
% Options::
% 'radian'     display angle in radians (default is degrees)
% 'fmt', f     use format string f for all numbers, (default %g)
% 'label',l    display the text before the transform
%
% Examples::
%        >> trprint2(T2)
%        t = (0,0), theta = -122.704 deg
%
%
% See also TRPRINT.

%## 2d homogeneous 


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

function out = trprint2(T, varargin)
    
    if ischar(T)
        % command form: trprint T
        trprint( evalin('caller', T) );
        return;
    end

    opt.fmt = [];
    opt.radian = false;
    opt.label = '';

    [opt,args] = tb_optparse(opt, varargin);

    s = '';

    if size(T,3) == 1
        if isempty(opt.fmt)
            opt.fmt = '%.3g';
        end
        s = tr2s(T, opt, args{:});
    else
        if isempty(opt.fmt)
            opt.fmt = '%8.2g';
        end        
        
        for i=1:size(T,3)
            % for each 4x4 transform in a possible 3D matrix
            s = char(s, tr2s(T(:,:,i), opt) );
        end
    end

    % if no output provided then display it
    if nargout == 0
        disp(s);
    else
        out = s;
    end
end

function s = tr2s(T, opt)
    % print the translational part if it exists
    if ~isempty(opt.label)
        s = sprintf('%8s: ', opt.label);
    else
        s = '';
    end
    if ~isrot2(T)
        s = strcat(s, sprintf('t = (%s),', vec2s(opt.fmt, transl2(T)')));
    end

    % print the angular part
    ang = atan2(T(2,1), T(1,1));
    if opt.radian
        s = strcat(s, ...
            sprintf(' %s rad', vec2s(opt.fmt, ang)) );
    else
        s = strcat(s, ...
            sprintf(' %s deg', vec2s(opt.fmt, ang*180.0/pi)) );
    end
end

function s = vec2s(fmt, v)
    s = '';
    for i=1:length(v)
        if abs(v(i)) < 1000*eps
            v(i) = 0;
        end
        s = [s, sprintf(fmt, v(i))];
        if i ~= length(v)
            s = [s, ', ']; % don't use strcat, removes trailing spaces
        end
    end
end
