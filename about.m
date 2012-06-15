%ABOUT Compact display of variable type
%
% ABOUT(X) displays a compact line that describes the class and dimensions of
% X.
%
% ABOUT X  as above but this is the command rather than functional form
%
% See also WHOS.
function about(var)
    
    if isstr(var)
        % invoked without parentheses
        w = evalin('caller', sprintf('whos(''%s'')', var));
        varname = var;
    else
        w = whos('var');
        varname = inputname(1);
    end
    
    if isempty(w)
        error(['cant find variable ' var])
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

