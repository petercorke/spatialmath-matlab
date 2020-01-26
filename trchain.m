%TRCHAIN Compound SE(3) transforms from string
%
% T = TRCHAIN(S) is a homogeneous transform (4x4) that results from
% compounding a number of elementary transformations defined by the string
% S.  The string S comprises a number of tokens of the form X(ARG) where
% X is one of Tx, Ty, Tz, Rx, Ry, or Rz.  ARG is an arbitrary MATLAB expression
% that can include constants or workspace variables. For example:
%
%          trchain('Tx(1) Rx(90) Ry(45) Tz(2)')
%
% is equivalent to computing
%
%        transl(1,0,0) * trotx(90, 'deg') * troty(45, 'deg') * transl(0,0,2)
%
% T = TRCHAIN(S, Q) as above but the expression for ARG can also contain
% a variable 'qJ' which selects the Jth value from the passed vector Q (1xN).
% For example:
%        trchain('Rx(q1)Tx(a1)Ry(q2)Ty(a3)Rz(q3)', [1 2 3])
%
% [T,TOK] = TRCHAIN(S ...) as above but return an array of tokens which can
% be passed in, instead of the string.
%
% T = TRCHAIN(TOK ...) as above but chain is defined by array of tokens
% instead of a string.
%
% Options::
% - 'deg'      all angular variables are in degrees (default radians)
% - 'qvar',V   treat the string V as the joint variable name rather than 'q'
%
% Notes::
% - Variables used in the string must exist in the caller workspace.
% - The string can contain arbitrary characters between the elements, for
%   example space, +, *, . or even |.
% - Works for symbolic variables in the workspace and/or passed in via the 
%   vector Q.
% - For symbolic operations that involve use of the value pi, make sure you
%   define it first in the workspace: pi = sym('pi');
% - The tokens are simply a parsed version of the input string and provide
%   some efficiency for repeated calls on the same chain.
%
% See also trchain2, trotx, troty, trotz, transl, SerialLink.trchain, ETS.

%## 3d homogeneous 

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


function [T,tokens_out] = trchain(s, varargin)
    
    opt.deg = false;
    opt.qvar = 'q';
    
    [opt,args] = tb_optparse(opt, varargin);
    
    if opt.deg
        unit = {'deg'};
    else
        unit = {};
    end
    
    if length(args) > 0
        q = args{1};
        if isa(q, 'symfun')
            q = formula(q);
        end
    else
        q = [];
    end
    
    if ischar(s) || isstring(s)
        % tokenise the input
        %    tok.op = {Tx, Ty, Tz, Rx, Ry, Rz}
        %    tok.arg = expression 
        %    tok.index = J if qJ is present, else is 0
        
        % get list of tokens
        tokens = regexp(char(s), '(?<op>R.?|T.)\((?<arg>[^\)]+)\)', 'names');
        
        % look for presence of a joint variable
        re = sprintf('%s(?<idx>[1-9][0-9]*)', opt.qvar);
        for i=1:length(tokens)
            qtok = regexp(tokens(i).arg, re, 'tokens');
            assert(length(qtok) <= 1, 'RTB:trchain:badarg', 'Only one qN variable allowed in any factor');
            if isempty(qtok)
                tokens(i).index = 0;
            else
                tokens(i).index = str2num(qtok{1}{1});
            end
        end
    else
        % token array provided, check it's good
        tokens = s;
        assert(all(isfield(s, {'arg', 'op', 'index'})), 'RTB:trchain:badarg', 'tokens are not a valid structure array');
    end
    
    T = eye(4,4);
    for token = tokens
        % for each token, evaluate the argument
        
        if token.index > 0
            % expression contains a joint variable
            % take value from the passed in vector q
            assert(token.index <= length(q), 'RTB:trchain:badarg', 'vector q has insufficient values');
            
            % substitute the joint value into the expression string
            if isa(q, 'sym')
                arg = subs(str2sym(token.arg), [opt.qvar num2str(token.index)], q(token.index));
            else
                expr = strrep(token.arg, [opt.qvar num2str(token.index)], num2str(q(token.index), '%g'));
                
                % evaluate and catch any errors
                try
                    arg = evalin('caller', expr);
                catch
                    error('RTB:trchain:badarg', 'can''t evaluate %s', expr);
                end
            end

        else
            % expression contains no joint variable
            % evaluate and catch any errors
            try
                arg = evalin('caller', token.arg);
            catch
                error('RTB:trchain:badarg', 'can''t evaluate %s', token.arg);
            end
        end
        
        % now compute the particular transform and update the transform chain
        switch token.op
            case 'Rx'
                T = T * trotx(arg, unit{:});
            case 'Ry'
                T = T * troty(arg, unit{:});
            case 'Rz'
                T = T * trotz(arg, unit{:});
            case 'Tx'
                T = T * transl(arg, 0, 0);
            case 'Ty'
                T = T * transl(0, arg, 0);
            case 'Tz'
                T = T * transl(0, 0, arg);
            otherwise
                error('RTB:trchain:badarg', 'unknown operator %s', token.op);
        end
    end
    
    if isa(q, 'sym')
        T = formula(T);
    end
    
    if nargout > 1
        tokens_out = tokens;
    end
end
