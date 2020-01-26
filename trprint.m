%TRPRINT Compact display of SE(3) homogeneous transformation
%
% TRPRINT(T, OPTIONS) displays the homogoneous transform (4x4) in a compact 
% single-line format.  If T is a homogeneous transform sequence then each 
% element is printed on a separate line.
%
% TRPRINT(R, OPTIONS) as above but displays the SO(3) rotation matrix (3x3).
%
% S = TRPRINT(T, OPTIONS) as above but returns the string.
%
% TRPRINT T OPTIONS is the command line form of above.

%
% Options::
% 'rpy'        display with rotation in ZYX roll/pitch/yaw angles (default)
% 'xyz'        change RPY angle sequence to XYZ
% 'yxz'        change RPY angle sequence to YXZ
% 'euler'      display with rotation in ZYZ Euler angles
% 'angvec'     display with rotation in angle/vector format
% 'radian'     display angle in radians (default is degrees)
% 'fmt', f     use format string f for all numbers, (default %g)
% 'label',l    display the text before the transform
% 'fid',F      send text to the file with file identifier F
%
% Examples::
%        >> trprint(T2)
%        t = (0,0,0), RPY/zyx = (-122.704,65.4084,-8.11266) deg
%
%        >> trprint(T1, 'label', 'A')
%               A:t = (0,0,0), RPY/zyx = (-0,0,-0) deg
%
%       >> trprint B euler
%       t = (0.629, 0.812, -0.746), EUL = (125, 71.5, 85.7) deg
%
% Notes::
% - If the 'rpy' option is selected, then the particular angle sequence can be
%   specified with the options 'xyz' or 'yxz' which are passed through to TR2RPY.
%  'zyx' is the default.
%
% See also TR2EUL, TR2RPY, TR2ANGVEC.

%## 3d homogeneous 

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

function out = trprint(T, varargin)
    
    opt.fmt = [];
    opt.mode = {'rpy', 'euler', 'angvec'};
    opt.radian = false;
    opt.label = [];
    opt.fid = 1;

    [opt,args] = tb_optparse(opt, varargin);
    
    if ischar(T)
        % command form: trprint T args
        opt.label = T;
        trprint( double(evalin('caller', T)), 'setopt', opt );
        return;
    end
    
    if isempty(opt.label)
        opt.label = inputname(1);
    end

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
            s = char(s, tr2s(T(:,:,i), opt, args{:}) );
        end
    end

    % if no output provided then display it
    if nargout == 0
        for row=s'
            fprintf(opt.fid, '%s\n', row);
        end
    else
        out = s;
    end
end

function s = tr2s(T, opt, varargin)
    % print the translational part if it exists
    if ~isempty(opt.label)
        s = sprintf('%8s: ', opt.label);
    else
        s = '';
    end
    if ~isrot(T)
        s = [s, sprintf('t = (%s),', vec2s(opt.fmt, transl(T)'))];
    end

    % print the angular part in various representations
    switch (opt.mode)
        case {'rpy', 'euler'}
            % angle as a 3-vector
            if strcmp(opt.mode, 'rpy')
                [ang,order] = tr2rpy(T, varargin{:});
                label = ['RPY/' order];
            else
                ang = tr2eul(T);
                label = 'EUL';
            end
            if opt.radian
                s = strcat(s, ...
                    sprintf(' %s = (%s) rad', label, vec2s(opt.fmt, ang)) );
            else
                s = strcat(s, ...
                    sprintf(' %s = (%s) deg', label, vec2s(opt.fmt, ang*180.0/pi)) );
            end
        case 'angvec'
            % as a vector and angle
            [th,v] = tr2angvec(T);
            if th == 0
                s = strcat(s, sprintf(' R = nil') );
            elseif opt.radian
                s = strcat(s, sprintf(' R = (%srad | %s)', ...
                    sprintf(opt.fmt, th), vec2s(opt.fmt, v)) );
            else
                s = strcat(s, sprintf(' R = (%sdeg | %s)', ...
                    sprintf(opt.fmt, th*180.0/pi), vec2s(opt.fmt,v)) );
            end
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
