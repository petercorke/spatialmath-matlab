%TRCHAIN Compound SE(3) transforms from string
%
% T = TRCHAIN(S, Q) is a homogeneous transform (4x4) that results from
% compounding a number of elementary transformations defined by the string
% S.  The string S comprises a number of tokens of the form X(ARG) where
% X is one of Tx, Ty, Tz, Rx, Ry, or Rz.  ARG is the name of a variable in
% MATLAB workspace or 'qJ' where J is an integer in the range 1 to N that
% selects the variable from the Jth column of the vector Q (1xN).
%
% For example:
%        trchain('Rx(q1)Tx(a1)Ry(q2)Ty(a3)Rz(q3)', [1 2 3])
%
% is equivalent to computing:
%        trotx(1) * transl(a1,0,0) * troty(2) * transl(0,a3,0) * trotz(3)
%
% Notes::
% - Variables list in the string must exist in the caller workspace.
% - The string can contain spaces between elements, or on either side of ARG.
% - Works for symbolic variables in the workspace and/or passed in via the 
%   vector Q.
% - For symbolic operations that involve use of the value pi, make sure you
%   define it first in the workspace: pi = sym('pi');
%
%
% See also trchain2, trotx, troty, trotz, transl, SerialLink.trchain, ETS.

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


function T = trchain(s, q)
    
    if nargin == 1
        q = [];
    end
    
    if isa(q, 'symfun')
        q = formula(q);
    end
    % s = 'Rx(q1)Tx(a1)Ry(q2)Tx(a3)Rz(q3)Tx(a3)';
    
    tokens = regexp(char(s), '\s*(?<op>R.?|T.)\(\s*(?<arg>[A-Za-z0-9\-][A-Za-z0-9+\-\*/]*)\s*\)\s*', 'names');
    
    T = eye(4,4);
    joint = 1;
    
    for token = tokens
        
        % get the argument for this transform element
        if token.arg(1) == 'q'
            % from the passed in vector q
            
            try
                arg = q(joint);
            catch
                error('SMTB:trchain:badarg', 'vector q has insufficient values');
            end
            joint = joint+1;
        else            % or the workspace
            
            try
                arg = evalin('caller', token.arg);
            catch
                error('SMTB:trchain:badarg', 'variable %s does not exist', token.arg);
            end
        end
        
        % now evaluate the element and update the transform chain
        
        switch token.op
            case 'Rx'
                T = T * trotx(arg);
            case 'Ry'
                T = T * troty(arg);
            case 'Rz'
                T = T * trotz(arg);
            case 'Tx'
                T = T * transl(arg, 0, 0);
            case 'Ty'
                T = T * transl(0, arg, 0);
            case 'Tz'
                T = T * transl(0, 0, arg);
            otherwise
                error('SMTB:trchain:badarg', 'unknown operator %s', token.op);
        end
    end
    
    if isa(q, 'symfun')
        T = formula(T);
    end
