
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

function smtbcheck
    % display current versions of MATLAB
    year = version('-release');
    fprintf('You are using:\n - MATLAB release %s\n', year);
    
    % check how old it is
    today = datevec(now);
    age =  today(1) - str2num(year(1:4));
    if age >= 2
        fprintf('     ** this is at least %d years old, you may have issues\n', age);
    end
    
    %% display versions of toolboxes (use unique RTB and MVTB functions)
    
    %SMTB
    smtb = checktoolbox('RTBPose', 'SMTB');
    
    
    %% check for shadowed files
    
    shadows = {};
    
    shadows = checkshadow(shadows, 'rotx');
    shadows = checkshadow(shadows, 'roty');
    shadows = checkshadow(shadows, 'rotz');
    shadows = checkshadow(shadows, 'angdiff');
    
    if ~isempty(shadows)
        fprintf('\n*** Some Toolbox files are "shadowed" and will cause problems with the use of this toolbox ***\n');
        for s = shadows'
            fprintf(' - toolbox function "%s" is shadowed by the MathWorks file %s\n', s{1}, s{2} );
        end
        fprintf('Use the pathtool function to:\n 1. move the Toolbox containing the shadowed function to the top of the path\n 2. Delete the problematic MathWorks toolbox (do you really need it?)\n 3. Move the problematic MathWorks toolbox to the end of the path\n')
    end
end

function present = checktoolbox(funcname, tbname)
    
    present = true;
    
    p = which(funcname);
    if isempty(p)
        % toolbox not present
        present = false;
        return
    end
    
    % display versions of toolbox
    fprintf('%s (%s)\n', short2long(tbname), p);
    v = ver( p );
    
    where = [];
    if strcmp(tbname, 'SMTB')
        % special case for SMTB, is it standalone or part of RTB or MVTB?
        if contains(p, 'vision', 'IgnoreCase', true)
            where = 'included with MVTB';
        elseif contains(p, 'robot', 'IgnoreCase', true)
            where = 'included with RTB';
        end
    end
    if isempty(where)
        if findstr(p, 'Add-Ons')
            where = 'Add-Ons (installed from mltbx file or Add-On Explorer)';
        else
            if exist( fullfile(p, '.git'), 'dir' )
                where = 'local (git clone)';
            else
                where = 'local (zip install)';
            end
        end
    end
    
    if ~isempty(v)
        fprintf(' - %s%s, %s\n', tbname, v.Version, v.Date);
    end
    fprintf(' - %s\n', where);
end


function out = checkshadow(shadows, funcname)
    % add to the list if function lives below MATLAB root
    funcpath = which(funcname); % find first instance in path
    out = shadows;
    if startsWith(funcpath, matlabroot) || startsWith(funcpath, 'built-in')
        out = [out; {funcname, which(funcname)}];
    end
end

function fullname = short2long(shortname)
    switch shortname
        case 'RTB'
            fullname = 'Robotics Toolbox for MATLAB';
        case 'MVTB'
            fullname = 'Machine Vision Toolbox for MATLAB';
        case 'SMTB'
            fullname = 'Spatial Math Toolbox for MATLAB';
    end
end

function p = getpath(funcname)
    funcpath = which(funcname);
    k = strfind( funcpath, [filesep funcname]);
    p = funcpath(1:k-1);
end