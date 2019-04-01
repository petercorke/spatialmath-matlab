%SE2 Representation of 2D rigid-body motion
%
% This subclasss of RTBPose is an object that represents rigid-body motion in 2D. 
% Internally this is a 3x3 homogeneous transformation matrix (3x3) belonging to 
% the group SE(2).
%
% Constructor methods::
%  SE2          general constructor
%  SE2.exp      exponentiate an se(2) matrix  
%  SE2.rand     random transformation
%  new          new SE2 object
%
% Display and print methods::
%  animate          ^graphically animate coordinate frame for pose
%  display          ^print the pose in human readable matrix form
%  plot             ^graphically display coordinate frame for pose
%  print            ^print the pose in single line format
%
% Group operations::
%  *            ^mtimes: multiplication (group operator, transform point)
%  /            ^mrdivide: multiply by inverse
%  ^            ^mpower: exponentiate (integer only): 
%  inv          inverse
%  prod         ^product of elements
%
% Methods::
%  det          determinant of matrix component
%  eig          eigenvalues of matrix component
%  log          logarithm of rotation matrix
%  inv          inverse
%  simplify*    apply symbolic simplication to all elements
%  interp       interpolate between poses
%  theta        rotation angle
%
% Information and test methods::
%  dim          ^returns 2
%  isSE         ^returns true
%  issym        ^test if rotation matrix has symbolic elements
%  SE2.isa      test if matrix is SE(2)
%
% Conversion methods::
%  char*         convert to human readable matrix as a string
%  SE2.convert   convert SE2 object or SE(2) matrix to SE2 object
%  double        convert to rotation matrix
%  R             convert to rotation matrix
%  SE3           convert to SE3 object with zero translation
%  SO2           convert rotational part to SO2 object
%  T             convert to homogeneous transformation matrix
%  Twist         convert to Twist object
%  t             get.t: convert to translation column vector
%
% Compatibility methods::
%  isrot2       ^returns false
%  ishomog2     ^returns true
%  tr2rt        ^convert to rotation matrix and translation vector
%  t2r          ^convert to rotation matrix
%  transl2      ^translation as a row vector  
%  trprint2     ^print single line representation
%  trplot2      ^plot coordinate frame
%  tranimate2   ^animate coordinate frame
%
% ^ inherited from RTBPose class.
%
% See also SO2, SE3, RTBPose.

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


classdef SE2 < SO2
    
    properties (Dependent = true)
        t
    end
    
    methods
        
        
        function obj = SE2(varargin)
            %SE2.SE2 Construct an SE(2) object
            %
            % Constructs an SE(2) pose object that contains a 3x3 homogeneous transformation
            % matrix.
            %
            % T = SE2() is the identity element, a null motion.
            %
            % T = SE2(X, Y) is an object representing pure translation defined by X and Y.
            %
            % T = SE2(XY) is an object representing pure translation defined by XY
            % (2x1). If XY (Nx2) returns an array of SE2 objects, corresponding to
            % the rows of XY.
            %
            % T = SE2(X, Y, THETA) is an object representing translation, X and Y, and
            % rotation, angle THETA.
            %
            % T = SE2(XY, THETA) is an object representing translation, XY (2x1), and
            % rotation, angle THETA.
            %
            % T = SE2(XYT) is an object representing translation, XYT(1) and XYT(2),
            % and rotation angle XYT(3). If XYT (Nx3) returns an array of SE2 objects, corresponding to
            % the rows of XYT.
            %
            % T = SE2(T) is an object representing translation and rotation defined by
            % the SE(2) homogeneous transformation matrix T (3x3).  If T (3x3xN) returns an 
            % array (1xN) of SE2 objects, corresponding to the third index of T.
            %
            % T = SE2(R) is an object representing pure rotation defined by the
            % SO(2) rotation matrix R (2x2)
            %
            % T = SE2(R, XY) is an object representing rotation defined by the
            % orthonormal rotation matrix R (2x2) and position given by XY (2x1)
            %
            % T = SE2(T) is a copy of the SE2 object T. If T (Nx1) returns an array of SE2 objects,
            % corresponding to the index of T.
            %
            % Options::
            % 'deg'         Angle is specified in degrees
            %
            % Notes::
            % - Arguments can be symbolic
            % - The form SE2(XY) is ambiguous with SE2(R) if XY has 2 rows, the second form is assumed.
            % - The form SE2(XYT) is ambiguous with SE2(T) if XYT has 3 rows, the second form is assumed.
            % - R and T are checked to be valid SO(2) or SE(2) matrices.
            
            opt.deg = false;
            
            [opt,args] = tb_optparse(opt, varargin);
            
            if opt.deg
                scale = pi/180.0;
            else
                scale = 1;
            end
            
            % if any of the arguments is symbolic the result will be symbolic
            if any( cellfun(@(x) isa(x, 'sym'), args) )
                obj.data = sym(obj.data);
            end
            
            obj.data = eye(3,3);

            switch length(args)
                case 0
                    % null motion
                    return
                case 1
                    % 1 argument
                    a = args{1};
                    
                    if isvec(a, 2)
                        % (t)
                        obj.data = [ 1 0 a(1); 0 1 a(2); 0 0 1];
                        
                    elseif isvec(a, 3)
                        % ([x y th])
                        a = a(:);
                        obj.data(1:2,1:2) = rot2(a(3)*scale);
                        obj.t = a(1:2);
                        
                    elseif SO2.isa(a)
                        % (R)
                        obj.data(1:2,1:2) = a;
                        
                    elseif SE2.isa(a)
                        % (T)
                        for i=1:size(a, 3)
                            obj(i).data = a(:,:,i);
                        end
                    elseif isa(a, 'SE2')
                        % (SE2)
                        for i=1:length(a)
                            obj(i).data = a(i).data;
                        end
                        
                    elseif any( numcols(a) == [2 3] )
                        for i=1:numrows(a)
                            obj(i) = SE2(a(i,:));
                        end
                        return
                    else
                        error('SMTB:SE2:badarg', 'unknown arguments');
                    end
                    
                case 2
                    % 2 arguments
                    a = args{1}; b = args{2};
                    if isscalar(a) && isscalar(b)
                        % (x,y)
                        obj.data = [ 1 0 a; 0 1 b; 0 0 1];
                    elseif isvec(a,2) && isscalar(b)
                        % ([x y], th)
                        obj.data = [ rot2(b*scale) a(:); 0 0 1];
                    elseif SO2.isa(a) && isvec(b,2)
                        % (R, t)
                        obj.data = [a b(:); 0 0 1];
                    else
                        error('SMTB:SE3:badarg', 'unknown arguments');
                    end
                    
                case 3
                    % 3 arguments
                    a = args{1}; b = args{2}; c = args{3};
                    if isscalar(a) && isscalar(b) && isscalar(c)
                        % (x, y, th)
                        obj.data = [ rot2(c*scale) [a b]'; 0 0 1];
                    else
                        error('SMTB:SE3:badarg', 'unknown arguments');
                    end
                otherwise
                    error('SMTB:SE3:badarg', 'unknown arguments');
                    
            end
            
            % add the last row if required
%             if numrows(obj.data) == 2
%                 obj.data = [obj.data; 0 0 1];
%             end
            assert(all(size(obj(1).data) == [3 3]), 'SMTB:SE2:SE2', 'created wrong size data element');
            %% HACK
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  GET AND SET
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function t = get.t(obj)
            %SE2.t  Get translational component
            %
            % P.t is a column vector (2x1) representing the translational component of
            % the rigid-body motion described by the SE2 object P.
            %
            % Notes::
            % - If P is a vector the result is a MATLAB comma separated list, in this
            %   case use P.transl().
            %
            % See also SE2.transl.
            t = obj.data(1:2,3);
        end
        
        function o = set.t(obj, t)
            %SE2.t  Set translational component
            %
            % P.t = TV sets the translational component of the rigid-body motion
            % described by the SE2 object P to TV (2x1).
            %
            % Notes::
            % - TV can be a row or column vector.
            % - If TV contains a symbolic value then the entire matrix becomes
            %   symbolic.
            
            if isa(t, 'sym') && ~isa(obj.data, 'sym')
                obj.data = sym(obj.data);
            end
            obj.data(1:2,3) = t;
            o = obj;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  conversion methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        function v = xyt(obj)
            %SE2.xyt  Extract configuration
            %
            % XYT = P.xyt() is a column vector (3x1) comprising the minimum three
            % configuration parameters of this rigid-body motion: translation (x,y)
            % and rotation theta.

            
            % TODO VECTORISE
            v = obj.t;
            v(3) = atan2(obj.data(2,1), obj.data(1,1));
        end        
        
        function [tx,ty] = transl(obj)
            %SE2.t  Get translational component
            %
            % TV = P.transl() is a row vector (1x2) representing the translational component of
            % the rigid-body motion described by the SE2 object P.  If P is a vector of
            % objects (1xN) then TV (Nx2) will have one row per object element.
            
            if nargout == 1 || nargout == 0
                tx = [obj.t]';
            else
                t = obj.t;
                tx = t(1);
                ty = t(2);
            end
        end
        
        function T = T(obj)
            %SE2.T  Get homogeneous transformation matrix
            %
            % T = P.T() is the homogeneous transformation matrix (3x3) associated with the
            % SE2 object P, and has zero translational component.  If P is a vector
            % (1xN) then T (3x3xN) is a stack of homogeneous transformation matrices, with the third
            % dimension corresponding to the index of P.
            %
            % See also SO2.T.
            for i=1:length(obj)
                T(:,:,i) = obj(i).data;
            end
        end
        
        function t = SE3(obj)
            %SE2.SE3 Lift to 3D
            %
            % Q = P.SE3() is an SE3 object formed by lifting the rigid-body motion
            % described by the SE2 object P from 2D to 3D.  The rotation is about the
            % z-axis, and the translation is within the xy-plane.
            %
            % See also SE3.
            t = SE3();
            t.data(1:2,1:2) = obj.data(1:2,1:2);
            t.data(1:2,4) = obj.data(1:2,3);
        end
        
        function out = SO2(obj)
            %SE2.SO2  Extract SO(2) rotation
            %
            % Q = SO2(P) is an SO2 object that represents the rotational component of
            % the SE2 rigid-body motion.
            %
            % See also SE2.R.
            
            out = SO2( obj.R );
        end
        
        function tw = Twist(obj)
            %SE2.Twist  Convert to Twist object
            %
            % TW = P.Twist() is the equivalent Twist object.  The elements of the twist are the unique
            % elements of the Lie algebra of the SE2 object P.
            %
            % See also SE2.log, Twist.
            tw = Twist( obj.log );
        end
        
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  SE(2) OPERATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function it = inv(obj)
            %SE2.inv  Inverse of SE2 object
            %
            % Q = inv(P) is the inverse of the SE2 object P.
            %
            % Notes::
            %  - This is formed explicitly, no matrix inverse required.
            %  - This is a group operator: input and output in the SE(2) group. 
            %  - P*Q will be the identity group element (zero motion, identity matrix).
 
            it = SE2( obj.R', -obj.R'*obj.t);
        end

        function S = log(obj)
            %SE2.log  Lie algebra
            %
            % se2 = P.log() is the Lie algebra corresponding to the SE2 object P. It is
            % an augmented skew-symmetric matrix (3x3).
            %
            % See also SE2.Twist, logm, skewa, vexa.
            S = logm(obj.data);
        end
        
        function Ti = interp(obj1, varargin)
            %SE2.interp Interpolate between SO2 objects
            %
            % P1.interp(P2, s) is an SE2 object which is an interpolation
            % between poses represented by SE2 objects P1 and P2.  s varies from 0
            % (P1) to 1 (P2). If s is a vector (1xN) then the result will be a vector
            % of SE2 objects.
            %
            % Notes::
            % - It is an error if S is outside the interval 0 to 1.
            %
            % See also SO2.angle.
            
            if isa(varargin{1}, 'SE2')
                % interp(SE2, SE2, s)  interpolate between given values
                obj2 = varargin{1};
                varargin = varargin(2:end);
                try
                    Ti = SE2( trinterp2(obj1.T, obj2.T, varargin{:}) );
                catch me
                    switch me.identifier
                        case 'SMTB:trinterp2:badarg'
                            throw( MException('SMTB:SE2:interp:badarg', 'value of S outside interval [0,1]') );
                        otherwise
                            rethrow(me);
                    end
                end
            else
                % interp(SE2, s)  interpolate between null and given value
                try
                    Ti = SE2( trinterp2( obj1.T, varargin{:}) );
                catch me
                    switch me.identifier
                        case 'SMTB:trinterp2:badarg'
                            throw( MException('SMTB:SE2:interp:badarg', 'value of S outside interval [0,1]') );
                        otherwise
                            rethrow(me);
                    end
                end
            end
            
        end

  
        function print(obj, varargin)
            for T=obj
                theta = atan2(T.data(2,1), T.data(1,1)) * 180/pi;
                fprintf('t = (%.4g, %.4g), theta = %.4g deg\n', T.t, theta);
            end
        end
        
        function n = new(obj, varargin)
            %SE2.new  Construct a new object of the same type
            %
            % P2 = P.new(X) creates a new object of the same type as P, by invoking the SE2 constructor on the matrix
            % X (3x3).
            %
            % P2 = P.new() as above but defines a null motion.
            %
            % Notes::
            %  - Serves as a dynamic constructor.
            %  - This method is polymorphic across all RTBPose derived classes, and
            %    allows easy creation of a new object of the same class as an existing
            %    one without needing to explicitly determine its type.
            %
            % See also SE3.new, SO3.new, SO2.new.
            
            n = SE2(varargin{:});
        end
        
    end
    
    methods (Static)
        % Static factory methods for constructors from exotic representations
        
        function obj = exp(s)
            %SE2.exp  Construct SE2 from Lie algebra
            %
            % SE2.exp(SIGMA) is the SE2 rigid-body motion corresponding to the se(2) 
            % Lie algebra element SIGMA (3x3).
            %
            % SE3.exp(TW) as above but the Lie algebra is represented
            % as a twist vector TW (1x1).
            %
            % Notes::
            %  - TW is the non-zero elements of X.
            %
            % Reference::
            % - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p25-31.
            %
            %
            % See also trexp2, skewa.

            obj = SE2( trexp2(s) );
        end
        

        
        function T = convert(tr)
            %SE2.check  Convert to SE2
            %
            % Q = SE2.convert(X) is an SE2 object equivalent to X where X is either
            % an SE2 object, or an SE(2) homogeneous transformation matrix (3x3).
            if isa(tr, 'SE2')
                T = tr;
            elseif SE2.isa(tr)
                T = SE2(tr);
            else
                error('expecting an SE2 or 3x3 matrix');
            end
        end
        
        function h = isa(tr, rtest)
            %SE2.ISA Test if matrix is SE(2)
            %
            % SE2.isa(T) is true (1) if the argument T is of dimension 3x3 or 3x3xN, else
            % false (0).
            %
            % SE2.isa(T, true) as above, but also checks the validity of the rotation
            % sub-matrix.
            %
            % Notes::
            %  - This is a class method.
            %  - The first form is a fast, but incomplete, test for a transform in SE(3).
            %  - There is ambiguity in the dimensions of SE2 and SO3 in matrix form.
            %
            % See also SO3.ISA, SE2.ISA, SO2.ISA, ishomog2.
            
            d = size(tr);
            if ndims(tr) >= 2
                h =  all(d(1:2) == [3 3]);
                
                if h && nargin > 1 && ~isa(tr, 'sym')
                    h = SO3.isa( tr(1:2,1:2) );
                    h = h && all(tr(4,:) == [0 0 0 1]);  % test the bottom row
                end
            else
                h = false;
            end
        end
                
        function T = rand()
            %SE2.rand Construct a random SE(2) object
            %
            % SE2.rand() is an SE2 object with a uniform random translation and a
            % uniform random orientation.  Random numbers are in the interval [-1 1] 
            % and rotations in the interval [-pi pi].
            %
            % See also RAND.
            T = SE2(rand(1,3)*diag([2 2 2*pi]) + [-1 -1 -pi]);
        end
    end
end
