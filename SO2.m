%SO2 Representation of 2D rotation
%
% This subclasss of RTBPose is an object that represents rotation in 2D.
% Internally this is a 2x2 orthonormal matrix belonging to the group SO(2).
%
% Constructor methods::
%  SO2          general constructor
%  SO2.exp      exponentiate an so(2) matrix
%  SO2.rand     random orientation
%  new          new SO2 object from instance
%
% Display and print methods::
%  animate     ^graphically animate coordinate frame for pose
%  display     ^print the pose in human readable matrix form
%  plot        ^graphically display coordinate frame for pose
%  print       ^print the pose in single line format
%
% Group operations::
%  *            ^mtimes: multiplication (group operator, transform point)
%  /            ^mrdivide: multiply by inverse
%  ^            ^mpower: exponentiate (integer only)
%  inv          ^inverse rotation
%  prod         ^product of elements
%
% Methods::
%  det          determinant of matrix value (is 1)
%  eig          ^eigenvalues of matrix value
%  interp       interpolate between rotations
%  log          logarithm of rotation matrix
%  simplify     ^apply symbolic simplication to all elements
%  subs         ^symbolic substitution
%  vpa          ^symbolic variable precision arithmetic
%
% Information and test methods::
%  dim          ^returns 2
%  isSE         ^returns false
%  issym        ^test if rotation matrix has symbolic elements
%  SO2.isa      test if matrix is SO(2)
%
% Conversion methods::
%  char          ^convert to human readable matrix as a string
%  SO2.convert   convert SO2 object or SO(2) matrix to SO2 object
%  double        ^convert to rotation matrix
%  theta         rotation angle
%  R             convert to rotation matrix
%  SE2           convert to SE2 object with zero translation
%  T             convert to homogeneous transformation matrix with zero translation
%
% Compatibility methods::
%  ishomog2     ^returns false
%  isrot2       ^returns true
%  tranimate2   ^animate coordinate frame
%  trplot2      ^plot coordinate frame
%  trprint2     ^print single line representation
%
% Operators::
%  +           ^plus: elementwise addition, result is a matrix
%  -           ^minus: elementwise subtraction, result is a matrix
%  ==          ^eq: test equality
%  ~=          ^ne: test inequality
%
% ^ inherited from RTBPose class.
%
% See also SE2, SO3, SE3, RTBPose.

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

classdef SO2 < RTBPose
    
    
    properties (Dependent = true)
        %R
    end
    
    methods
        
        function obj = SO2(t, deg)
            %SO2.SO2  Construct SO2 object
            %
            % P = SO2() is the identity element, a null rotation.
            %
            % P = SO2(THETA) is an SO2 object representing rotation of THETA radians.
            % If THETA is a vector (N) then P is a vector of objects, corresponding to
            % the elements of THETA.
            %
            % P = SO2(THETA, 'deg') as above but with THETA degrees.
            %
            % P = SO2(R) is an SO2 object formed from the rotation 
            % matrix R (2x2).
            %
            % P = SO2(T) is an SO2 object formed from the rotational part 
            % of the homogeneous transformation matrix T (3x3).
            %
            % P = SO2(Q) is an SO2 object that is a copy of the SO2 object Q.
            %
            % Notes::
            %  - For matrix arguments R or T the rotation submatrix is checked for validity.
            %
            % See also rot2, SE2, SO3.
            
            if nargin == 0
                % null rotation
                obj.data = eye(2,2);
            elseif isa(t, 'SO2')
                % copy an existing SO2 object
                obj.data = t.data;
            elseif isvector(t)
                % for specified angle
                if isfloat(t)
                    assert(isreal(t), 'argument cannot be complex');
                end
                t = t(:).';
                if nargin > 1 && strcmp(deg, 'deg')
                    t = t *pi/180;
                end
                for i=1:length(t)
                    th = t(i);
                    obj(i).data = [
                        cos(th)  -sin(th)
                        sin(th)   cos(th)
                        ];
                end

            elseif SO2.isa(t)
                % from a 2x2 matrix
                for i=1:size(t, 3)
                    x = t(:,:,i);
                    assert(SO2.isa(x, 'valid'), 'SMTB:SO2.SO2:badarg', 'matrix is not in SO(2)');
                    obj(i).data = x;
                end
            elseif SE2.isa(t)
                % from a 3x3 matrix
                for i=1:size(t, 3)
                    x = t(1:2,1:2,i);
                    assert(SO2.isa(x, 'valid'), 'SMTB:SO2.SO2:badarg', 'rotation submatrix is not in SO(2)');
                    obj(i).data = x;
                end  
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  GET AND SET
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        function RR = R(obj)
            %SO2.R  Get rotation matrix
            %
            % R = P.R() is the rotation matrix (2x2) associated with the SO2 object P.
            % If P is a vector (1xN) then R (2x2xN) is a stack of rotation matrices, with
            % the third dimension corresponding to the index of P.
            %
            % See also SO2.T.
            if ~issym(obj)
                RR = zeros(2,2,length(obj)); % prealloc so long as not symbolic
            end
            for i=1:length(obj)
                RR(:,:,i) = obj(i).data(1:2,1:2);
            end
        end
        
        function TT = T(obj)
            %SO2.T  Get homogeneous transformation matrix
            %
            % T = P.T() is the homogeneous transformation matrix (3x3) associated with the
            % SO2 object P, and has zero translational component.  If P is a vector
            % (1xN) then T (3x3xN) is a stack of rotation matrices, with the third
            % dimension corresponding to the index of P.
            %
            % See also SO2.T.
            TT = zeros(3,3,length(obj));
            for i=1:length(obj)
                TT(1:2,1:2,i) = obj(i).data(1:2,1:2);
                TT(3,3,i) = 1;
            end
        end
        
        function th = angle(obj)
            %SO2.theta  Rotation angle
            %
            % P.angle() is the rotation angle, in radians [-pi,pi), associated with the
            % SO2 object P.
            %
            % See also atan2.
            th = atan2(obj.data(2,1), obj.data(1,1));
        end
  
        function th = theta(obj)
            %SO2.theta  Rotation angle
            %
            % P.theta() is the rotation angle, in radians, associated with the
            % SO2 object P.
            %
            % Notes::
            %  - Deprecated, use angle() instead.
            %
            % See also SO2.angle.
            
            % TODO: remove
            
            warning('SMTB:SO2:deprecfun', 'the theta() method is deprecated use angle()');
            th = obj.angle();
        end

                
        function s = char(obj)
            %SO2.char Convert to string
            %
            % P.char() is a string containing rotation matrix elements.
            %
            % See also RTB.display.            
            s = num2str(obj.data, 4);
        end
        
        function T = SE2(obj)
            %SO2.SE2 Convert to SE2 object
            %
            % P.SE2() is an SE2 object formed from the rotational component of the
            % SO2 object P and with a zero translational component.
            %
            % See also SE2.
            T = SE2( r2t(obj.data) );
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  OPERATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ir = inv(obj)
            %SO2.inv  Inverse
            %
            % Q = inv(P) is an SO2 object representing the inverse of the SO2 object P.  
            %
            % Notes::
            %  - This is a group operator: input and output in the SO(2) group.
            %  - This is simply the transpose of the underlying matrix.
            %  - P*Q will be the identity group element (zero rotation, identity matrix).
            ir = SO2(obj.data');
        end
        
        function d = det(obj)
            %SO2.inv  Determinant
            %
            % det(P) is the determinant of the SO2 object P and should always be +1.
            d = det(obj.R);
        end
        
        function varargout = eig(obj, varargin)
            %SO2.eig  Eigenvalues and eigenvectors
            %
            % E = eig(P) is a column vector containing the eigenvalues of the
            % underlying rotation matrix.
            %
            % [V,D] = eig(P) produces a diagonal matrix D of eigenvalues and 
            % a full matrix V whose columns are the corresponding eigenvectors  
            % such that A*V = V*D.
            %
            % See also eig.
            [varargout{1:nargout}] = eig(obj.data, varargin{:});
        end
            
        function S = log(obj)
            %SO2.log  Logarithm
            %
            % so2 = P.log() is the Lie algebra corresponding to the SO2 object P. It is
            % a skew-symmetric matrix (2x2).
            %
            % See also SO2.exp, Twist, logm, vex, skew.

            S = logm(obj.data);
        end
        
        %         function o = set.R(obj, data)
        %             if isa(data, 'sym')
        %                 obj.data = data;
        %             else
        %             obj.data(1:2,1:2) = data;
        %             end
        %             o = obj;
        %         end

        
        function R = interp(obj1, obj2, s)
            %SO2.interp Interpolate between rotations
            %
            % P1.interp(P2, s) is an SO2 object representing interpolation
            % between rotations represented by SO2 objects P1 and P2.  s varies from 0
            % (P1) to 1 (P2). If s is a vector (1xN) then the result will be a vector
            % of SO2 objects.
            %
            % P1.interp(P2,N) as above but returns a vector (1xN) of SO2 objects
            % interpolated between P1 and P2 in N steps.
            %
            % Notes::
            %  - It is an error if any element of S is outside the interval 0 to 1.
            %
            % See also SO2.angle.
            assert(all(s>=0 & s<=1), 'SMTB:SO2:interp:badarg', 's must be in the interval [0,1]');
            
            th1 = obj1.angle; th2 = obj2.angle;
            
            R = SO2( th1 + s*(th2-th1) );
        end
        
        
        function R = trinterp2(obj, varargin)
            R = obj.interp(varargin{:});
        end
        

        
        function n = new(obj, varargin)
            %SO2.new  Construct a new object of the same type
            %
            % Create a new object of the same type as the RTBPose derived instance object.
            %
            % P.new(X) creates a new object of the same type as P, by invoking the SO2 constructor on the matrix
            % X (2x2).
            %
            % P.new() as above but assumes an identity matrix.
            %
            % Notes::
            %  - Serves as a dynamic constructor.
            %  - This method is polymorphic across all RTBPose derived classes, and
            %   allows easy creation of a new object of the same class as an existing
            %   one without needing to explicitly determine its type.
            %
            % See also SE3.new, SO3.new, SE2.new.
            
            n = SO2(varargin{:});
        end

        function r = times(p, q)
            r = mtimes(p, q);
        end

    end
    methods (Static)
        % Static factory methods for constructors from exotic representations
        
        function obj = exp(s)
            %SO2.exp  Construct SO2 from Lie algebra
            %
            % R = SO3.exp(X) is the SO2 rotation corresponding to the so(2) 
            % Lie algebra element SIGMA (2x2).
            %
            % R = SO3.exp(TW) as above but the Lie algebra is represented
            % as a twist vector TW (1x1).
            %
            % Notes::
            %  - TW is the non-zero elements of X.
            %
            % Reference::
            %  - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p25-31.
            %
            % See also trexp2, skewa.
            obj = SO2( trexp2(s) );
        end
       
        
        function R = convert(tr)
            %SO2.convert  Convert value to SO2
            %
            % Q = SO2.convert(X) is an SO2 object equivalent to X where X is either 
            % an SO2 object, an SO(2) rotation matrix (2x2), an SE2 object, or an
            % SE(2) homogeneous transformation matrix (3x3).
            if isa(tr, 'SO2')
                R = SO2(tr);        % is SO2 or SE2, enforce it being an SO2
            elseif SO2.isa(tr)
                R = SO2(tr);        % is 2x2
            elseif SE2.isa(tr)
                R = SO2( t2r(tr) );  % is 3x3, assume SE(2), take rotational part
            else
                error('SMTB:SO2:convert:badarg', 'expecting an SO2, 2x2, SE3 or 3x3');
            end
        end
        
        
        function h = isa(r, dtest)
            %SO2.ISA Test if matrix belongs to SO(2)
            %
            % SO2.ISA(T) is true (1) if the argument T is of dimension 2x2 or 2x2xN, else
            % false (0).
            %
            % SO2.ISA(T, true) as above, but also checks the validity of the rotation
            % matrix, ie. that its determinant is +1.
            %
            % Notes::
            %  - The first form is a fast, but incomplete, test for a transform in SO(2).
            %
            % See also SO3.ISA, SE2.ISA, SE2.ISA, ishomog2.
            d = size(r);
            if (isfloat(r) || isa(r, 'sym') ) && ndims(r) >= 2
                h =  all(d(1:2) == [2 2]);
                
                if h && nargin > 1 && ~isa(r, 'sym')
                    h = abs(det(r) - 1) < eps;
                end
            else
                h = false;
            end
        end
        
        function T = rand()
            %SO2.rand Construct a random SO(2) object
            %
            % SO2.rand() is an SO2 object where the angle is drawn from a uniform 
            % random orientation. Random numbers are in the interval 0 to 2pi.
            %
            % See also RAND.
            T = SO2(rand(1,1)*2*pi);
        end
    end
end

