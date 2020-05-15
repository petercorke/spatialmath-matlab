%RTBPose Superclass for SO2, SO3, SE2, SE3
%
% This abstract class provides common methods for the 2D and 3D orientation and pose
% classes: SO2, SE2, SO3 and SE3.
%
% Display and print methods::
%  animate          graphically animate coordinate frame for pose
%  display          print the pose in human readable matrix form
%  plot             graphically display coordinate frame for pose
%  print            print the pose in single line format
%
% Group operations::
%  *           mtimes: multiplication within group, also transform vector
%  /           mrdivide: multiplication within group by inverse
%  prod        mower: product of elements
%
% Methods::
%
%  dim         dimension of the underlying matrix
%  isSE        true for SE2 and SE3
%  issym       true if value is symbolic
%  simplify    apply symbolic simplification to all elements
%  vpa         apply vpa to all elements
%
% % Conversion methods::
%  char        convert to human readable matrix as a string
%  double      convert to real rotation or homogeneous transformation matrix
%
% Operators::
%  +           plus: elementwise addition, result is a matrix
%  -           minus: elementwise subtraction, result is a matrix
%  ==          eq: test equality
%  ~=          ne: test inequality
%
% Compatibility methods::
% A number of compatibility methods give the same behaviour as the
% classic RTB functions:
%
%  tr2rt       convert to rotation matrix and translation vector
%  t2r         convert to rotation matrix
%  tranimate   animate coordinate frame
%  trprint     print single line representation
%  trprint2    print single line representation
%  trplot      plot coordinate frame
%  trplot2     plot coordinate frame
%
% Notes::
% - This is a handle class.
% - RTBPose subclasses can be used in vectors and arrays.
% - Multiplication and division with normalization operations are performed
%   in the subclasses.
% - SO3 is polymorphic with UnitQuaternion making it easy to change
%   rotational representations.
%
% See also SO2, SO3, SE2, SE3.

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


classdef (Abstract) RTBPose
    
    properties(Access=protected, Hidden=true)
        data   % this is a 2x2, 3x3 or 4x4 matrix, possibly symbolic
    end
    
    properties
        ref    % string, name of reference coordinate frame
        target % string, name of target coordinate frame
    end
    
    
    methods
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  INFORMATION METHODS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function n = dim(obj)
            %RTBPose.dim  Dimension
            %
            % N = P.dim() is the dimension of the matrix representing the RTBPose
            % subclass instance P.  It is 2 for SO2, 3 for SE2 and SO3, and 4 for SE3.
            
            n = size(obj(1).data, 2);
        end
        
        function t = isSE(T)
            %RTBPose.isSE  Test if rigid-body motion
            %
            % P.isSE() is true if P is an instance of the RTBPose sublass SE2 or SE3.
            
            s = class(T);
            t = s(2) == 'E';
        end
        
        function t = issym(obj)
            %RTBPose.issym  Test if pose is symbolic
            %
            % P.issym() is true if the RTBPose subclass instance P has symbolic rather 
            % than real values.
            t = isa(obj(1).data, 'sym');
        end
        
%         function e = isequal(obj1, obj2)
%             %ISEQUAL Test quaternion element equality
%             %
%             % ISEQUAL(P1,P2) is true if the RTBPose subclass instances P1 
%%            %
%             % Notes::
%             % - Used by test suite verifyEqual in addition to eq().
%             % - Invokes eq().
%             %
%             % See also Quaternion.eq.
%             e = isa(obj2, classname(obj1)) && ...
%                 length(obj1) == length(obj2) &&
%                 alleq(q1, q2);
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  OPERATORS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function v = plus(a, b)
            %RTBPose.plus Add poses
            %
            % P1+P2 is the elementwise summation of the matrix elements of the
            % RTBPose subclass instances P1 and P2.  The result is a native matrix not 
            % the input class type since the result of addition is not in the group.
            
            assert(isa(b, class(a)) && length(a) == length(b), 'SMTB:RTBPose:plus:badarg', 'operands don''t conform');
            for i=1:length(a)
                v(:,:,i) = a(i).data + b(i).data;
            end
        end
        
        function v = minus(a, b)
            %RTBPose.minus Subtract poses
            %
            % P1-P2 is the elementwise difference of the matrix elements of the two
            % poses.  The result is a matrix not the input class type since the result
            % of subtraction is not in the group.
            
            assert(isa(b, class(a)) && length(a) == length(b), 'SMTB:RTBPose:minus:badarg', 'operands don''t conform');
            for i=1:length(a)
                v(:,:,i) = a(i).data - b(i).data;
            end
        end
        
%        function v = uminus(a)
%            %RTBPose.uminus Unary minus of poses
%            %
%            % -P is the elementwise negation of the matrix elements of the 
%            % pose.
%            
%            for i=1:length(a)
%                v = a.new(-a(i).data);
%            end
%        end
        
        function out = simplify(obj)
            %RTBPose.simplify Symbolic simplification
            %
            % P2 = P.simplify() applies symbolic simplification to each element of
            % internal matrix representation of the RTBPose subclass instance P.
            %
            % See also simplify.
            out = obj;
            if isa(obj(1).data, 'sym')
                for k=1:length(obj)
                    % simplify every element of data
                    for i=1:numel(obj.data)
                        out(k).data(i) = simplify( obj(k).data(i) );
                    end
                end
            end
        end
        
        function out = vpa(obj, D)
            %RTBPose.vpa Variable precision arithmetic
            %
            % P2 = P.vpa() numerically evaluates each element of
            % internal matrix representation of the RTBPose subclass instance P.
            %
            % P2 = P.vpa(D) as above but with D decimal digit accuracy.
            %
            % Notes::
            %  - Values of symbolic variables are taken from the workspace.
            %
            % See also vpa, simplify.
            out = obj;
            if nargin == 1
                D = digits;
            end
            if isa(obj(1).data, 'sym')
                for k=1:length(obj)
                    % simplify every element of data
                    for i=1:numel(obj.data)
                        out(k).data(i) = vpa( obj(k).data(i), D );
                    end
                end
            end
        end
        
        function e = isequal(obj1, obj2)
            e = eq(obj1, obj2);
        end
        
        function e = eq(obj1, obj2)
            e = false;
            
            if ~isa(obj2, class(obj2)) || ~(length(obj1) == length(obj2))
                return;
            end
            
            assert(length(obj1) == length(obj2), 'SMTB:RTBPose:eq', 'arrays must be same size');
            e = zeros(size(obj1), 'logical');
            for i=1:length(obj1)
                e(i) = all(all(abs([obj1(i).data] - [obj2(i).data]) < 10*eps));
            end
            
        end
        
        function e = ne(obj1, obj2)
            e = ~eq(obj1, obj2);
        end
        

        function e = mpower(obj1, n)
            %RTBPose.mpower Exponential of pose
            %
            % P^N is an RTBPose subclass instance equal to RTBPose subclass instance P raised
            % to the integer power N.  It is equivalent of compounding P with itself N-1 times.
            %
            % Notes::
            % - N can be 0 in which case the result is the identity element.
            % - N can be negative which is equivalent to the inverse of P^(-N).
            %
            % See also RTBPose.power, RTBPose.mtimes, RTBPose.times.
            assert(isscalar(n) && isreal(n) && floor(n) == n, 'SMTB:Pose', 'exponent must be a real integer');
            e = obj1.new( double(obj1)^n);
        end
        
        function e = power(obj1, n)
            %RTBPose.power Exponential of pose
            %
            % P.^N is the exponential of P where N is an integer, followed by normalization.  It is equivalent of compounding
            % the rigid-body motion of P with itself N-1 times.
            %
            % Notes::
            % - N can be 0 in which case the result is the identity matrix.
            % - N can be negative which is equivalent to the inverse of P.^abs(N).
            %
            % See also RTBPose.mpower, RTBPose.mtimes, RTBPose.times.
            assert(isscalar(n) && isreal(n) && floor(n) == n, 'SMTB:Pose', 'exponent must be a real integer');
            x = double(obj1)^n;
            switch obj1.dim
                case 2
                    e = obj1.new(x);
                case 3
                    e = obj1.new( trnorm(x) );
            end
        end
        
        function e = transpose(obj1, obj2)
            error('SMTB:Pose', 'transpose operator not supported by RTBPose subclass object')
        end
        function e = ctranspose(obj1, obj2)
            error('SMTB:Pose', 'transpose operator not supported by RTBPose subclass object')
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  COMPOSITION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = mtimes(obj, a)
            %RTBPose.mtimes  Compound pose objects
            %
            % R = P*Q is an RTBPose subclass instance representing the composition of the
            % RTBPose subclass instance P by the RTBPose subclass instance Q.
            %
            % If either, or both, of P or Q are vectors, then the result is a vector.
            %  - if P is a vector (1xN) then R is a vector (1xN) such that R(i) = P(i)*Q.
            %  - if P is a vector (1xN) then R is a vector (1xN) such that R(i) = P*Q(i).
            %  - if both P and Q are vectors (1xN) then R is a vector (1xN) such 
            %    that R(i) = P(i)*Q(i).
            %
            % W = P*V is a column vector (2x1) which is the transformation of the
            % column vector V (2x1) by the matrix representation of the RTBPose
            % subclass instance P.
            %
            % P can be a vector and/or V can be a matrix, a columnwise set of vectors:
            %  - if P is a vector (1xN) then W is a matrix (2xN) such that W(:,i) = P(i)*V.
            %  - if V is a matrix (2xN) V is a matrix (2xN) then W is a matrix (2xN) such
            %    that W(:,i) = P*V(:,i).
            %  - if P is a vector (1xN) and V is a matrix (2xN) then W is a matrix (2xN)
            %    such that W(:,i) = P(i)*V(:,i).
            %
            % Notes::
            % - Computed by matrix multiplication of their equivalent matrices.
            %
            % See also RTBPose.mrdivide.

            if strcmp(class(obj), class(a))
                % obj * obj
                obj1 = obj(1);
                out = repmat(obj1, 1, max(length(obj),length(a)));
                if length(obj) == length(a)
                    % do objvector*objvector and objscalar*objscalar case
                    for i=1:length(obj)
                        out(i) = obj1.new( obj(i).data * a(i).data);
                    end
                elseif length(obj) == 1
                    % objscalar*objvector case
                    for i=1:length(a)
                        out(i) = obj1.new( obj.data * a(i).data);
                    end
                elseif length(a) == 1
                    % objvector*objscalar case
                    for i=1:length(obj)
                        out(i) = obj1.new( obj(i).data * a.data);
                    end
                else
                    error('SMTB:RTBPose:badops', 'invalid operand lengths to * operator');
                end
                
            elseif isa(obj, 'RTBPose') && isnumeric(a)
                % obj * vectors (nxN), result is nxN
                assert(isreal(a), 'SMTB:RTBPose:*', 'matrix must be real');
                q2 = double(a); % force to double
                
                obj1 = obj(1);
                n = numrows(obj1.data);

                if obj1.isSE
                    % is SE(n) convert to homogeneous form
                    assert(numrows(a) == n-1, 'SMTB:RTBPose:badops', 'LHS should be matrix with %d rows', n-1);
                    a = [a; ones(1, numcols(a))];
                else
                    assert(numrows(a) == n, 'SMTB:RTBPose:badops', 'LHS should be matrix with %d rows', n);
                end
                
                out = zeros(n, max(length(obj),numcols(a)));  % preallocate space
                
                if length(obj) == numcols(a)
                    % do objvector*vector and objscalar*scalar case
                    for i=1:length(obj)
                        out(:,i) = obj(i).data * a(:,i);
                    end
                elseif length(obj) == 1
                    % objscalar*vector case
                    for i=1:length(obj)
                        out = obj.data * a;
                    end
                elseif numcols(a) == 1
                    % objvector*scalar case
                    for i=1:length(obj)
                        out(:,i) = obj(i).data * a;
                    end
                else
                    error('SMTB:RTBPose:badops', 'unequal vector lengths to * operator');
                end
                
                if obj1.isSE
                    % is SE(n) convert to homogeneous form
                    out = out(1:end-1,:);
                end
            elseif isa(obj, 'SE2') && isa(a, 'polyshape')
                    % special case, planar rigid body transform of a polyshape
                    out = polyshape( (obj * a.Vertices')' );
            elseif isa(obj, 'SE3') && isa(a, 'Plucker')
                A = [obj.R -skew(obj.t); zeros(3,3) obj.R];
                out = A * a;   % invokes mtimes Plucker.mtimes
            else
                error('SMTB:RTBPose:badops', 'operands to * cannot be composed');
            end
        end
        
        function out = prod(obj)
            %RTBPose.prod Compound array of poses
            %
            % P.prod() is an RTBPose subclass instance representing the product (composition) of the
            % successive elements of P (1xN).
            %
            % Note::
            % - Composition is performed with the .* operator, ie. the product is 
            %   renormalized at every step.
            %
            % See also RTBPose.times.
            out = obj(1);
            
            for i=2:length(obj)
                out = out .* obj(i);
            end
        end
        
        function out = mrdivide(obj, a)
            %RTBPose.mrdivide  Compound SO2 object with inverse
            %
            % R = P/Q is an RTBPose subclass instance representing the composition of the
            % RTBPose subclass instance P by the inverse of the RTBPose subclass instance Q.
            %
            % If either, or both, of P or Q are vectors, then the result is a vector.
            %  - if P is a vector (1xN) then R is a vector (1xN) such that R(i) = P(i)/Q.
            %  - if P is a vector (1xN) then R is a vector (1xN) such that R(i) = P/Q(i).
            %  - if both P and Q are vectors (1xN) then R is a vector (1xN) such 
            %    that R(i) = P(i)/Q(i).
            %
            % Notes::
            % - Computed by matrix multiplication of their equivalent matrices with 
            %   the second one inverted.
            %
            % See also RTBPose.mtimes.
            
            obj1 = obj(1);
            n = obj1.dim;
            
            if strcmp(class(obj1), class(a))
                % obj / obj
                out = repmat(obj1, 1, max(length(obj),length(a)));
                if length(obj) == length(a)
                    % do vector/vector and scalar/scalar case
                    for i=1:length(obj)
                        out(i) = obj1.new( obj(i).data * inv(a(i).data));
                    end
                elseif length(obj) == 1
                    % scalar/vector case
                    for i=1:length(a)
                        out(i) = obj1.new( obj.data * inv(a(i).data) );
                    end
                elseif length(a) == 1
                    % vector/scalar case
                    for i=1:length(obj)
                        out(i) = obj1.new( obj(i).data * inv(a.data));
                    end
                else
                    error('SMTB:RTBPose:badops', 'unequal vector lengths to / operator');
                end 
            else
                error('SMTB:RTBPose:badops', 'invalid operand types to / operator');
            end
        end
  
        function v = subs(obj, old, new)
            %RTBPose.subs Symbolic substitution
            %
            % T = subs(T, old, new) replaces old with new in the symbolic
            % transformation T.
            %
            % See also: subs

            v = obj.new();  % clone the input
            v.data = subs(v.data, old, new);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  COMPATABILITY/CONVERSION METHODS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [R,t] = tr2rt(obj)
            %tr2rt  Split rotational and translational components  (compatibility)
            %
            % [R,t] = tr2rt(P) is the rotation matrix and translation vector
            % corresponding to the SE2 or SE3 instance P.
            %
            % See also tr2rt.
            n = numcols(obj.data);
            
            assert(isSE(obj), 'only applicable to SE2/3 class');
            if length(obj) > 1
                R = zeros(3,3,length(obj));
                t = zeros(length(obj), 3);
                for i=1:length(obj)
                    R(:,:,i) = obj(i).R;
                    t(i,:) = obj(i).t';
                end
            else
                R = obj.R;
                t = obj.t;
            end
        end
        
        function R = t2r(obj)
            %t2r  Get rotation matrix  (compatibility)
            %
            % t2r(P) is a native matrix corresponding to the rotational
            % component of the SE2 or SE3 instance P.
            %
            % See also t2r.

            n = numcols(obj.data);
            
            assert(isSE(obj), 'only applicable to SE2/3 class');
            if length(obj) > 1
                R = zeros(3,3,length(obj));
                for i=1:length(obj)
                    R(:,:,i) = obj(i).R;
                end
            else
                R = obj.R;
            end
        end
        
                
        function d = double(obj)
            %RTBPose.double  Convert to matrix
            %
            % T = P.double() is a native matrix representation of the RTBPose 
            % subclass instance P, either a rotation matrix or a homogeneous 
            % transformation matrix.
            %
            % If P is a vector (1xN) then T will be a 3-dimensional array (MxMxN).
            %
            % Notes::
            % - If the pose is symbolic the result will be a symbolic matrix.
            
            if ~isa(obj(1).data, 'sym')
                d = zeros( [size(obj(1).data) length(obj)] );
            end
            for i=1:length(obj)
                d(:,:,i) = obj(i).data;
            end
        end
        
        
        function out = trprint(obj, varargin)
            %TRPRINT Compact display of 3D rotation or transform (compatibility)
            %
            % trprint(P, OPTIONS) displays the RTBPose subclass instance P in a 
            % compact single-line format.  If P is a vector then each element is 
            % printed on a separate line.
            %
            % Notes::
            %  - see trprint for details of options.
            %  - P can be instances of SO3 or SE3.
            %
            % See also RTBPose.print, trprint.
            args = {varargin{:}, 'label', inputname(1)};
            
            if nargout == 0
                print(obj, args{:});
            else
                out = print(obj, args{:});
            end
        end
        
        function out = trprint2(obj, varargin)
            %TRPRINT2 Compact display of 2D rotation or transform (compatibility)
            %
            % trprint2(P, OPTIONS) displays the RTBPose subclass instance P in a 
            % compact single-line format.  If P is a vector then each element is 
            % printed on a separate line.
            %
            % Notes::
            %  - see trprint for details of options.
            %  - P can be instances of SO2 or SE2.
            %
            % See also RTBPose.print, trprint2.
            if nargout == 0
                print(obj, varargin{:});
            else
                out = print(obj, varargin{:});
            end  
        end
        
        function varargout = trplot(obj, varargin)
            %TRPLOT Draw a 3D coordinate frame (compatibility)
            %
            % trplot(P, OPTIONS) draws a 3D coordinate frame represented by RTBPose
            % subclass instance P.
            %
            % Notes::
            %  - see trplot for details of options.
            %  - P can be instances of SO3 or SE3.
            %
            % See also RTBPose.plot, trplot.
            [varargout{1:nargout}] = obj.plot(varargin{:});
        end
        
        function varargout = trplot2(obj, varargin)
            %TRPLOT2 Draw a 2D coordinate frame (compatibility)
            %
            % trplot2(P, OPTIONS) draws a 2D coordinate frame represented by RTBPose
            % subclass instance P.
            %
            % Notes::
            %  - see trplot for details of options.
            %  - P can be instances of SO2 or SE2.
            %
            % See also RTBPose.plot, trplot2.
            [varargout{1:nargout}] = obj.plot(varargin{:});
        end
        
        function tranimate(obj, varargin)
            %TRANIMATE Animate a 3D coordinate frame (compatibility)
            %
            % TRANIMATE(P1, P2, OPTIONS) animates a 3D coordinate frame moving between
            % RTBPose subclass instances P1 and pose P2.
            %
            % TRANIMATE(P, OPTIONS) animates a 2D coordinate frame moving from the identity pose
            % to the RTBPose subclass instance P.
            %
            % TRANIMATE(PV, OPTIONS) animates a trajectory, where PV is a vector of
            % RTBPose subclass instances.
            %
            % Notes::
            %  - see tranimate for details of options.
            %  - P, P1, P2, PV can be instances of SO3 or SE3.
            %
            % See also RTBPose.animate, tranimate.
            
            obj.animate(varargin{:});
        end
        
        function tranimate2(obj, varargin)
            %TRANIMATE2 Animate a 2D coordinate frame (compatibility)
            %
            % TRANIMATE2(P1, P2, OPTIONS) animates a 2D coordinate frame moving between
            % RTBPose subclass instances P1 and pose P2.
            %
            % TRANIMATE2(P, OPTIONS) animates a 2D coordinate frame moving from the identity pose
            % to the RTBPose subclass instance P.
            %
            % TRANIMATE2(PV, OPTIONS) animates a trajectory, where PV is a vector of
            % RTBPose subclass instances.
            %
            % Notes::
            %  - see tranimate2 for details of options.
            %  - P, P1, P2, PV can be instances of SO2 or SE2.
            %
            % See also RTBPose.animate, tranimate.
            obj.animate(varargin{:});
        end
        
        function v = isrot(obj)
            %ISROT Test if SO3 class (compatibility)
            %
            % ISROT(R) is true (1) if R is of class SO3.
            %
            % See also ISROT.
            v = obj.dim == 3 && ~obj.isSE;
        end
        
        function v = isrot2(obj)
            %ISROT2 Test if SO2 class (compatibility)
            %
            % ISROT2(R) is true (1) if R is of class SO2.
            %
            % See also ISROT2.

            v = obj.dim == 2 && ~obj.isSE;
        end
        
        function v = ishomog(obj)
            %ISHOMOG Test if SE3 class (compatibility)
            %
            % ISHOMOG(T) is true (1) if T is of class SE3.
            %
            % See also ISHOMOG.
            v = obj.dim == 4 && obj.isSE;
        end
        
        function v = ishomog2(obj)
            %ISHOMOG2 Test if SE2 class (compatibility)
            %
            % ISHOMOG2(T) is true (1) if T is of class SE2.
            %
            % See also ISHOMOG2.
            v = obj.dim == 3 && obj.isSE;
        end

        function v = isvec(obj, n)
            %ISVEC Test if vector (compatibility)
            %
            % ISVEC(T) is always false.
            %
            % See also ISVEC.
            v = false;
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  DISPLAY METHODS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function display(obj)
            %RTBPose.display Display pose in matrix form
            %
            % P.display() displays the matrix elements for the RTBPose instance P to
            % the console. If P is a vector (1xN) then matrices are displayed sequentially.
            %
            % Notes::
            % - This method is invoked implicitly at the command line when the result
            %   of an expression is an RTBPose subclass object and the command has no trailing
            %   semicolon.
            % - If the function cprintf is found is used to colorise the matrix:
            %   rotational elements in red, translational in blue. 
            %   See https://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
%
            %
            % See also SO2, SO3, SE2, SE3.
            
            try  % for Octave
                loose = strcmp( get(0, 'FormatSpacing'), 'loose');
            catch
                loose = 0;
            end
            if loose
                disp(' ');
            end
            
            obj.render(inputname(1));  % the hard work done in render
        end
        
        function disp(obj)
            disp( char(obj) );
        end
        
        function s2 = char(obj)
            %RTBPose.char Convert to string
            %
            % s = P.char() is a string showing RTBPose matrix elements as
            % a matrix.
            %
            % See also RTBPose.display.
            s = num2str(obj.data, '%10.4g'); %num2str(obj.data, 4);
            for i=1:numrows(s);
                s2(i,:) = ['    ', s(i,:)];
            end
        end

        
        function out = print(obj, varargin)
            %RTBPose.print Compact display of pose
            %
            % P.print(OPTIONS) displays the RTBPose subclass instance P in a compact
            % single-line format.  If P is a vector then each element is printed on 
            % a separate line.
            %
            % Example::
            %          T = SE3.rand()
            %          T.print('rpy', 'xyz')  % display using XYZ RPY angles
            %
            % Notes::
            % - Options are passed through to trprint or trprint2 depending on the object
            %   type.
            %
            % See also trprint, trprint2.
            if isa(obj, 'SO2') || isa(obj, 'SE2')
                printfn = @trprint2;
            else
                printfn = @trprint;
            end
            if nargout == 0
                for T=obj
                    printfn(double(T), varargin{:});
                end
            else
                out = '';
                
                for T=obj
                    out = char(out, printfn(double(T), varargin{:}));
                end
            end
        end
        
        function animate(obj, varargin)
            %RTBPose.animate Animate a coordinate frame
            %
            % RTBPose.animate(P1, P2, OPTIONS) animates a 3D coordinate frame moving from
            % RTBPose P1 to RTBPose P2.
            %
            % RTBPose.animate(P, OPTIONS) animates a coordinate frame moving from the identity pose
            % to the RTBPose P.
            %
            % RTBPose.animate(PV, OPTIONS) animates a trajectory, where PV is a vector of
            % RTBPose subclass objects.
            %            %
            % Options::
            %  'fps', fps    Number of frames per second to display (default 10)
            %  'nsteps', n   The number of steps along the path (default 50)
            %  'axis',A      Axis bounds [xmin, xmax, ymin, ymax, zmin, zmax]
            %  'movie',M     Save frames as files in the folder M
            %  'cleanup'     Remove the frame at end of animation
            %  'noxyz'       Don't label the axes
            %  'rgb'         Color the axes in the order x=red, y=green, z=blue
            %  'retain'      Retain frames, don't animate
            % 
            % Additional options are passed through to tranimate or tranimate2.
            %
            % See also tranimate, tranimate2.
            
            % invoke classic functions
            if length(varargin) > 0 && isa(varargin{1}, 'RTBPose')
                % tranimate(T1, T2, args)
                switch class(obj)
                    case 'SO2'
                        tranimate2(obj.R, varargin{1}.R, varargin{2:end});
                        
                    case 'SE2'
                        tranimate2(obj.T, varargin{1}.T, varargin{2:end});
                        
                    case 'SO3'
                        tranimate(obj.R, varargin{1}.R, varargin{2:end});
                        
                    case 'SE3'
                        tranimate(obj.T, varargin{1}.T, varargin{2:end});
                end
            else
                % tranimate(T1, args)
                switch class(obj)
                    case 'SO2'
                        tranimate2(obj.R, varargin{:});
                        
                    case 'SE2'
                        tranimate2(obj.T, varargin{:});
                        
                    case 'SO3'
                        tranimate(obj.R, varargin{:});
                        
                    case 'SE3'
                        tranimate(obj.T, varargin{:});
                end
            end
            

        end
        
        function varargout = plot(obj, varargin)
            %TRPLOT Draw a coordinate frame (compatibility)
            %
            % trplot(P, OPTIONS) draws a 3D coordinate frame represented by P which is
            % SO2, SO3, SE2 or SE3.
            %
            % Compatible with matrix function trplot(T).
            %
            % Options are passed through to trplot or trplot2 depending on the object
            % type.
            %
            % See also trplot, trplot2.
             
             
            switch class(obj)
                case 'SO2'
                    [varargout{1:nargout}] = trplot2(obj.R, varargin{:});
                    
                case 'SE2'
                    [varargout{1:nargout}] = trplot2(obj.T, varargin{:});
                    
                case 'SO3'
                    [varargout{1:nargout}] = trplot(obj.R, varargin{:});
                    
                case 'SE3'
                    [varargout{1:nargout}] = trplot(obj.T, varargin{:});
            end
        end
        
        
    end
    
    methods (Access=private)
        function render(obj, varname)
            
            if isa(obj(1).data, 'sym')
                % use MATLAB default disp() function for symbolic object
                disp(obj.data);
            else
                % else render the elements with specified format and color
                fmtR = '%10.4f';
                fmtt = '%10.4g';
                fmt0 = '%10.0f';
                if exist('cprintf')
                    print = @(color, fmt, value) cprintf(color, fmt, value);
                else
                    print = @(color, fmt, value) fprintf(fmt, value);
                end
                switch class(obj)
                    case 'SO2',  nr = 2; nt = 0;
                    case 'SE2',  nr = 2; nt = 1;
                    case 'SO3',  nr = 3; nt = 0;
                    case 'SE3',  nr = 3; nt = 1;
                end
                
                for i=1:length(obj)
                    M = obj(i).data;
                    if length(obj) > 1
                        fprintf('\n%s(%d) = \n', varname, i);
                    else
                        fprintf('\n%s = \n', varname);
                        
                    end
                    M(abs(M)<1000*eps) = 0;
                    
                    for row=1:nr
                        for col=1:(nr+nt)
                            if col <= nr
                                % rotation matrix
                                v = M(row,col);
                                
                                if fix(v) == v
                                    print('Errors', fmt0, v); % red
                                else
                                    print('Errors', fmtR, v); % red
                                end
                            else
                                % translation
                                print('Keywords', fmtt, M(row,col)); % blue
                            end
                        end
                        fprintf('\n');
                    end
                    % last row
                    if nt > 0
                        for col=1:(nr+nt)
                            print('Text', fmt0, M(nr+nt,col));
                        end
                        fprintf('\n');
                    end
                end
            end
        end
    end
    
end
