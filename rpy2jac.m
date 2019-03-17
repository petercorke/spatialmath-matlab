%RPY2JAC Jacobian from RPY angle rates to angular velocity
%
% J = RPY2JAC(RPY, OPTIONS) is a Jacobian matrix (3x3) that maps ZYX roll-pitch-yaw angle 
% rates to angular velocity at the operating point RPY=[R,P,Y].
%
% J = RPY2JAC(R, P, Y, OPTIONS) as above but the roll-pitch-yaw angles are passed
% as separate arguments.
%
% Options::
% 'xyz'     Use XYZ roll-pitch-yaw angles
% 'yxz'     Use YXZ roll-pitch-yaw angles
%
% Notes::
% - Used in the creation of an analytical Jacobian.
% - Angles in radians, rates in radians/sec.
%
% Reference::
% - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p232-3.
%
% See also eul2jac, rpy2r, SerialLink.jacobe.

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

function J = rpy2jac(r, varargin)

    opt.order = {'zyx', 'xyz', 'yxz'};
    [opt,args] = tb_optparse(opt, varargin);
    
    
        % unpack the arguments
    if numcols(r) == 3
		p = r(:,2);
		y = r(:,3);
		r = r(:,1);
	elseif nargin >= 3
        p = args{1};
        y = args{2};
    else
        error('SMTB:rpy2jac:badarg', 'bad arguments')
    end
    
    
    switch opt.order
    case 'xyz'
        J = [	
        sin(p)          0       1  
        -cos(p)*sin(y)  cos(y)  0  
        cos(p)*cos(y)   sin(y)  0  
        ];
    
     case 'zyx'
        J = [ 
            cos(p)*cos(y), -sin(y), 0
            cos(p)*sin(y),  cos(y), 0
            -sin(p),       0, 1
            ];
    
    case 'yxz'
        J = [
            cos(p)*sin(y),  cos(y), 0
            -sin(p),       0, 1
            cos(p)*cos(y), -sin(y), 0
            ];
    end
    
%{
    syms r p y rd pd yd wx wy wz real
    syms rt(t) pt(t) yt(t) 

    order = 'yxz'

    R = rpy2r(r, p, y, order);
    Rt = rpy2r(rt, pt, yt, order);
    dRdt = diff(Rt, t);
    dRdt = subs(dRdt, {diff(rt(t),t), diff(pt(t),t), diff(yt(t),t),}, {rd,pd,yd});
    dRdt = subs(dRdt, {rt(t),pt(t),yt(t)}, {r,p,y});
    dRdt = formula(dRdt)   % convert symfun to an array

    w = vex(dRdt * R');
    w = simplify(w)

    clear A
    rpyd = [rd pd yd];

    for i=1:3
        for j=1:3
            C = coeffs(w(i), rpyd(j));
            if length(C) == 1
                A(i,j) = 0;
            else
            A(i,j) = C(2);
            end
        end
    end

    A
%}
