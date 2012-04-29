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
    fprintf('%s [%s] : %s (%d bytes)\n', ...
        varname, w.class, ss, w.bytes);