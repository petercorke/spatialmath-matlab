
% display current versions of MATLAB
year = version('-release');
fprintf('You are using:\n - MATLAB release %s\n', year);

% check how old it is
today = datevec(now);
age =  today(1) - str2num(year(1:4));
if age >= 2
fprintf('     ** this is at least %d years old, you may have issues\n', age);
end
    
% display versions of toolboxes (use unique RTB and MVTB functions)
a = ver( getpath('tr2rpy') );
rtb = ~isempty(a);
if rtb
    fprintf(' - %s %s %s\n', a.Name, a.Version, a.Date);
end
a = ver( getpath('idisp') );
mvtb = ~isempty(a);
if mvtb
    fprintf(' - %s %s %s\n', a.Name, a.Version, a.Date);
end

% check for shadowed files
k = 0;
if rtb
    k = k + checkpath('rotx');
    k = k + checkpath('roty');
    k = k + checkpath('rotz');
end
if mvtb
    k = k + checkpath('im2col');
    k = k + checkpath('col2im');
    k = k + checkpath('angdiff');
end

if k > 0
    fprintf('Some Toolbox files are shadowed.  Use path tool to move the toolbox to the top of the path\n')
end

function k = checkpath(funcname)
    
    funcpath = which(funcname); % find first instance in path
    k = 0;

    good = {'rvc', 'Robotics Toolbox for MATLAB', 'Machine Vision Toolbox for MATLAB'};
    if exist(funcname)
        if all( cellfun(@(x) isempty(strfind(funcpath, x)), good) )
            fprintf('** Toolbox function %s is shadowed by %s\n', funcname, which(funcname) );
            k = 1;
        end
    end
end

function p = getpath(funcname)
        funcpath = which(funcname);
        k = strfind( funcpath, [filesep funcname]);
        p = funcpath(1:k-1);
end