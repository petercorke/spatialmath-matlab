%Quaternion Quaternion class
%
% A quaternion is 4-element mathematical object comprising a scalar s, and
% a vector v which can be considered as a pair (s,v).  In the Toolbox it is
% denoted by q = s <<vx, vy, vz>>.
%
% A quaternion of unit length can be used to represent 3D orientation and
% is implemented by the subclass UnitQuaternion.
%
% Constructors::
%  Quaternion        general constructor
%  Quaternion.pure   pure quaternion
%
% Display and print methods::
%  display           print in human readable form
%
% Group operations::
%  *       quaternion (Hamilton) product or elementwise multiplication by scalar
%  /       multiply by inverse or elementwise division by scalar
%  ^       exponentiate (integer only)
%  +       elementwise sum of quaternion elements 
%  -       elementwise difference of quaternion elements
%  conj    conjugate
%  exp     exponential
%  log     logarithm
%  inv     inverse
%  prod    product of elements
%  unit    unitized quaternion
%
% Methods::
%  inner             inner product
%  isequal           test for non-equality
%  norm              norm, or length
%
% Conversion methods::
%  char              convert to string
%  double            quaternion elements as 4-vector
%  matrix            quaternion as a 4x4 matrix
%
% Overloaded operators::
%  ==      test for quaternion equality
%  ~=      test for quaternion inequality
%
% Properties (read only)::
%  s         real part
%  v         vector part
%
% Notes::
% - This is reference (handle) class object
% - Quaternion objects can be used in vectors and arrays.
%
% References::
% - Animating rotation with quaternion curves, K. Shoemake,
%   in Proceedings of ACM SIGGRAPH, (San Fran cisco), pp. 245-254, 1985.
% - On homogeneous transforms, quaternions, and computational efficiency,
%   J. Funda, R. Taylor, and R. Paul,
%   IEEE Transactions on Robotics and Automation, vol. 6, pp. 382-388, June 1990.
% - Quaternions for Computer Graphics, J. Vince, Springer 2011.
% - Robotics, Vision & Control: Second Edition, P. Corke, Springer 2016; p44-45.
%
% See also UnitQuaternion.

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


% TODO
%  constructor handles R, T trajectory and returns vector
%  .r, .t on a quaternion vector??

classdef Quaternion
    
    properties (SetAccess = protected, GetAccess=public)
        s       % scalar part
        v       % vector part
    end
    
    
    methods
        
        function q = Quaternion(s, v)
            %Quaternion.Quaternion Construct a quaternion object
            %
            % Q = Quaternion(S, V) is a Quaternion formed from the scalar S and vector
            % part V (1x3).
            %
            % Q = Quaternion([S V1 V2 V3]) is a Quaternion formed by specifying directly its 4 elements.
            %
            % Q = Quaternion() is a zero Quaternion, all its elements are zero.
            %
            % Notes::
            % - The constructor is not vectorized, it cannot create a vector of Quaternions.
            
            if nargin == 0
                q.v = [0,0,0];
                q.s = 0;
            elseif isa(s, 'Quaternion')
                q.s = s.s;
                q.v = s.v;
            elseif nargin == 2 && isscalar(s) && isvec(v,3)
                q.s = s;
                q.v = v(:).';
            elseif nargin == 1 && isvec(s,4)
                s = s(:).';
                q.s = s(1);
                q.v = s(2:4);
            elseif isrot(s) || isa(s, 'SO3')
                error('SMTB:Quaternion:badarg', 'regular Quaternion is not equivalent to a rotation, use UnitQuaternion instead')
            else
                error ('SMTB:Quaternion:badarg', 'bad argument to quaternion constructor');
            end
            
        end

        function qo = set.s(q, s)
            %Quaternion.set.s Set scalar component
            %
            % Q.s = S sets the scalar part of the Quaternion object to S.
            assert(isa(s, 'sym') || ( isreal(s) && isscalar(s) ), 'SMTB:Quaternion:badarg', 's must be real scalar');
            
            qo = q;
            qo.s = s;
        end
        
        function qo = set.v(q, v)
            %Quaternion.set.v Set vector component
            %
            % Q.v = V sets the vector part of the Quaternion object to V (1x3).
            qo = q;
            if isa(v, 'symfun')
                qo.v = v;
            else
                assert(isvec(v,3), 'SMTB:Quaternion:badarg', 'v must be a real 3-vector');
                
                qo.v = v(:).';
            end
        end
        
%         function s = get.s(q)
%             s = [q.s]';
%         end
%         
%         function v = get.v(q)
%             [q.v]
%             v = reshape([q.v]', 3, [])';
%         end
            
        function display(q)
            %Quaternion.display Display quaternion
            %
            % Q.display() displays a compact string representation of the Quaternion's value
            % as a 4-tuple.  If Q is a vector then S has one line per element.
            %
            % Notes::
            % - This method is invoked implicitly at the command line when the result
            %   of an expression is a Quaternion object and the command has no trailing
            %   semicolon.
            % - The vector part is displayed with double brackets << 1, 0, 0 >> to
            %   distinguish it from a UnitQuaternion which displays as < 1, 0, 0 >
            % - If Q is a vector of Quaternion objects the elements are displayed on
            %   consecutive lines.
            %
            % See also Quaternion.char.
            
            loose = strcmp( get(0, 'FormatSpacing'), 'loose');
            if loose
                disp(' ');
            end
            disp([inputname(1), ' = '])
            if loose
                disp(' ');
            end
            disp(char(q))
            if loose
                disp(' ');
            end
        end
        
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% QUATERNION FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        function c = conj(q)
            %Quaternion.conj Conjugate of a quaternion
            %
            % Q.conj() is a Quaternion object representing the conjugate of Q.
            %
            % Notes::
            % - Conjugatation is the negation of the vector component.
            %
            % See also Quaternion.inv.
            c = q.new(q.s, -q.v);
        end
        
        function qi = inv(q)
            %Quaternion.inv Invert a quaternion
            %
            % Q.inv() is a Quaternion object representing the inverse of Q.
            %
            % Notes::
            % - If Q is a vector then an equal length vector of Quaternion objects
            %   is computed representing the elementwise inverse of Q.
            %
            % See also Quaternion.conj.
            
            for i=1:length(q)
                n2 = sum( q(i).double.^2 );
                qi(i) = Quaternion([q(i).s -q(i).v]/ n2);
            end
        end
        
        function qu = unit(q)
            %Quaternion.unit Unitize a quaternion
            %
            % QU = Q.unit() is a Quaternion with a norm of 1.  If Q is a vector (1xN) then
            % QU is also a vector (1xN).
            %
            % Notes::
            % - This is Quaternion of unit norm, not a UnitQuaternion object.
            %
            % See also Quaternion.norm, UnitQuaternion.
            
            for i=1:length(q)
                qu(i) = Quaternion( q(i).double / norm(q(i)) );
            end
        end
       
        function n = norm(q)
            %Quaternion.norm Quaternion magnitude
            %
            % Q.norm(Q) is the scalar norm or magnitude of the Quaternion Q.
            %
            % Notes::
            % - This is the Euclidean norm of the Quaternion written as a 4-vector.
            % - A unit-quaternion has a norm of one and is represented by the
            %   UnitQuaternion class.
            %
            % See also Quaternion.inner, Quaternion.unit, UnitQuaternion.
            
            n = colnorm(double(q)')';
        end
        
        function m = matrix(q)
            %Quaternion.matrix Matrix representation of Quaternion
            %
            % Q.matrix() is a matrix (4x4) representation of the Quaternion Q.
            %
            % Quaternion, or Hamilton, multiplication can be implemented as a
            % matrix-vector product, where the column-vector is the elements of a
            % second quaternion:
            %
            %          matrix(Q1) * double(Q2)'
            %
            % Notes::
            % - This matrix is not unique, other matrices will serve the purpose for
            %   multiplication, see https://en.wikipedia.org/wiki/Quaternion#Matrix_representations
            % - The determinant of the matrix is the norm of the Quaternion to the fourth power. 
            %
            % See also Quaternion.double, Quaternion.mtimes.
            m = [q.s    -q.v(1) -q.v(2) -q.v(3)
                 q.v(1)  q.s    -q.v(3)  q.v(2)
                 q.v(2)  q.v(3)  q.s    -q.v(1)
                 q.v(3) -q.v(2)  q.v(1)  q.s];
        end
        
        function n = inner(q1, q2)
            %Quaternion.inner Quaternion inner product
            %
            % V = Q1.inner(Q2) is the inner (dot) product of two vectors (1x4),
            % comprising the elements of Q1 and Q2 respectively.
            %
            % Notes::
            % - Q1.inner(Q1) is the same as Q1.norm().
            %
            % See also Quaternion.norm.
            
            n = double(q1)*double(q2)';
        end
        
        function out = log(q)
            %Quaternion.log Logarithm of quaternion
            %
            % Q.log() is the logarithm of the Quaternion Q.
            %
            % See also Quaternion.exp.
            
            assert(norm(q.v) > 20*eps, 'SMTB:Quaternion:log:badarg', 'Can''t compute log of Quaternion with zero length vector component');
            out = Quaternion( log(norm(q)), unit(q.v) * acos(q.s/norm(q)) );
        end
        
        function out = exp(q)
            %Quaternion.log Exponential of quaternion
            %
            % Q.log() is the logarithm of the Quaternion Q.
            %
            % See also Quaternion.exp.
           
            assert(norm(q.v) > 20*eps, 'SMTB:Quaternion:exp:badarg', 'Can''t compute exp of Quaternion with zero length vector component');
            out = exp(q.s) * Quaternion( cos(norm(q.v)), unit(q.v)*sin(norm(q.v)) );
        end
        
        function out = prod(q)
            %Quaternion.prod Product of quaternions
            %
            % prod(Q) is the product of the elements of the vector of Quaternion objects Q.
            %
            % See also Quaternion.mtimes, RTBPose.prod.
            out = q(1);
            for qq = q(2:end)
                out = out * qq;
            end
        end
       
        function out = spositive(q)
            out = q;
            if q.s < 0
                out.s = -out.s;
                out.v = -out.v;
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ARITHMETIC OPERATORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        function qp = mtimes(q1, q2)
            %Quaternion.mtimes Multiply a quaternion object
            %
            % Q1*Q2   is a Quaternion formed by the Hamilton product of two Quaternions.
            % Q*S     is the element-wise multiplication of Quaternion elements by the scalar S.
            % S*Q     is the element-wise multiplication of Quaternion elements by the scalar S.
            %
            % Notes::
            % - Overloaded operator '*'.
            % - If either, or both, of Q1 or Q2 are vectors, then the result is a vector.
            %   - if Q1 is a vector (1xN) then R is a vector (1xN) such that R(i) = Q1(i)*Q2.
            %   - if Q2 is a vector (1xN) then R is a vector (1xN) such that R(i) = Q1*Q2(i).
            %   - if both Q1 and Q2 are vectors (1xN) then R is a vector (1xN) such 
            %     that R(i) = Q1(i)*Q2(i).
            %
            % See also Quaternion.mrdivide, Quaternion.mpower.
            
            if isa(q1, 'Quaternion') && isa(q2, 'Quaternion')
                %QQMUL  Multiply quaternion by quaternion
                %
                % QQ = qqmul(Q1, Q2) is the product of two quaternions.
                
                if isa(q1, 'UnitQuaternion') && isa(q2, 'UnitQuaternion')
                    new = @UnitQuaternion.new;
                    newclass = 'UnitQuaternion';
                else
                    new = @Quaternion.new;
                    newclass = 'Quaternion';
                end
                if all(size(q1) == size(q2))
                    for i=1:length(q1)
                        % decompose into scalar and vector components
                        s1 = q1(i).s;  v1 = q1(i).v;
                        s2 = q2(i).s;  v2 = q2(i).v;
                        
                        % form the product
                        qp(i) = new([s1*s2-v1*v2.' s1*v2+s2*v1+cross(v1,v2)]);
                    end
                elseif isscalar(q1)
                    s1 = q1.s;  v1 = q1.v;
                    
                    for i=1:length(q2)
                        % decompose into scalar and vector components
                        s2 = q2(i).s;  v2 = q2(i).v;
                        
                        % form the product
                        qp(i) = new([s1*s2-v1*v2.' s1*v2+s2*v1+cross(v1,v2)]);
                    end
                elseif isscalar(q2)
                    s2 = q2.s;  v2 = q2.v;

                    for i=1:length(q1)
                        % decompose into scalar and vector components
                        s1 = q1(i).s;  v1 = q1(i).v;
                        
                        % form the product
                        qp(i) = new([s1*s2-v1*v2.' s1*v2+s2*v1+cross(v1,v2)]);
                    end
                else
                    error('SMTB:Quaternion:badarg', '* operand length/size mismatch');
                end
                
            elseif isa(q1, 'Quaternion') && isa(q2, 'double')
                
                %QSMUL  Multiply quaternion
                %
                % Q = qsmul(Q, S) multiply quaternion by real scalar.
                %
                assert(isscalar(q2), 'SMTB:Quaternion:badarg', 'quaternion-double product: must be a scalar');
                for i=1:length(q1)
                    qp(i) = Quaternion( double(q1(i))*q2);
                end
                
                
            elseif isa(q1, 'double') && isa(q2, 'Quaternion')
                
                %QSMUL  Multiply quaternion
                %
                % Q = qsmul(Q, S) multiply quaternion by real scalar.
                %
                
                assert(isscalar(q1), 'SMTB:Quaternion:badarg', 'quaternion-double product: must be a scalar');
                
                for i=1:length(q2)
                    qp(i) = Quaternion( double(q2(i))*q1);
                end
            else
                error('SMTB:Quaternion:badarg', 'quaternion product: incorrect right hand operand');
            end
        end
       
        function qq = mrdivide(q1, q2)
            %Quaternion.mrdivide Quaternion quotient.
            %
            % R = Q1/Q2   is a Quaternion formed by Hamilton product of Q1 and inv(Q2).
            % R = Q/S     is the element-wise division of Quaternion elements by the scalar S.
            %
            % Notes::
            % - Overloaded operator '/'.
            % - If either, or both, of Q1 or Q2 are vectors, then the result is a vector.
            %   - if Q1 is a vector (1xN) then R is a vector (1xN) such that R(i) = Q1(i)./Q2.
            %   - if Q2 is a vector (1xN) then R is a vector (1xN) such that R(i) = Q1./Q2(i).
            %   - if both Q1 and Q2 are vectors (1xN) then R is a vector (1xN) such 
            %     that R(i) = Q1(i)./Q2(i).
            %
            % See also Quaternion.mtimes, Quaternion.mpower, Quaternion.plus, Quaternion.minus.
            
            if isa(q1, 'Quaternion') && isa(q2, 'Quaternion')
                %QQDIV  Divide quaternion by quaternion
                %
                % QQ = qqdiv(Q1, Q2) is the quotient of two quaternions.
                
                if length(q1) == length(q2)
                    for i=1:length(q1)
                        
                        % form the quotient
                        qq(i) = q1(i) * inv(q2(i));
                    end
                elseif isscalar(q1)
                    
                    for i=1:length(q2)
                      
                        % form the quotient
                        qq(i) = q1 * inv(q2(i));
                    end
                elseif isscalar(q2)
                    
                    for i=1:length(q1)
                        
                        % form the quotient
                        qq(i) = q1(i) * inv(q2);
                    end
                else
                    error('SMTB:Quaternion:badarg', '/ operand length mismatch');
                end
                
            elseif isa(q1, 'Quaternion') && isa(q2, 'double')
                
                %QSDIV  Divide quaternion by scalar
                %
                % Q = qsdiv(Q, S) divide quaternion by real scalar.
                %
                
                assert(isscalar(q2), 'SMTB:Quaternion:badarg', 'quaternion-double quotient: must be a scalar');
                for i=1:length(q1)
                    qq(i) = Quaternion( double(q1(i))/q2);
                end
            
            else
                error('SMTB:Quaternion:badarg', 'quaternion quotient: incorrect right hand operand');
            end
        end

                
        function qp = mpower(q, p)
            %Quaternion.mpower Raise quaternion to integer power
            %
            % Q^N is the Quaternion Q raised to the integer power N.
            %
            % Notes::
            % - Overloaded operator '^'.
            % - N must be an integer, computed by repeated multiplication.
            %
            % See also Quaternion.mtimes.
            
            % check that exponent is an integer
            assert(p - floor(p) == 0, 'SMTB:Quaternion:badarg', 'quaternion exponent must be integer');
            
            if p == 0
                qp = q.new([1 0 0 0]);
            else
                qp = q;
                
                % multiply by itself so many times
                for i = 2:abs(p)
                    qp = qp * q;
                end
                
                % if exponent was negative, invert it
                if p<0
                    qp = inv(qp);
                end
            end
        end
        
        
        function qp = plus(q1, q2)
            %PLUS Add quaternions
            %
            % Q1+Q2 is a Quaternion formed from the element-wise sum of Quaternion elements.
            %
            % Q1+V  is a Quaternion formed from the element-wise sum of Q1 and the
            % vector V (1x4).
            %
            % Notes::
            % - Overloaded operator '+'.
            % - Effectively V is promoted to a Quaternion.
            %
            % See also Quaternion.minus.
            
            if isa(q2, 'Quaternion')
                qp = Quaternion(double(q1) + double(q2));
            elseif isvec(q2, 4)
                qp = Quaternion(q1);
                q2 = q2(:)';
                qp.s = qp.s + q2(1);
                qp.v = qp.v + q2(2:4);
            end
        end
        
        function qp = minus(q1, q2)
            %Quaternion.minus Subtract quaternions
            %
            % Q1-Q2 is a Quaternion formed from the element-wise difference of Quaternion elements.
            %
            % Q1-V  is a Quaternion formed from the element-wise difference of Q1 and the
            % vector V (1x4).
            %
            % Notes::
            % - Overloaded operator '-'.
            % - Effectively V is promoted to a Quaternion.
            %
            % See also Quaternion.plus.
            
            if isa(q2, 'Quaternion')
                qp = Quaternion(double(q1) - double(q2));
            elseif isvec(q2, 4)
                qp = Quaternion(q1);
                q2 = q2(:)';
                qp.s = qp.s - q2(1);
                qp.v = qp.v - q2(2:4);
            end
        end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% RELATIONAL OPERATORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        function e = isequal(q1, q2)
            %ISEQUAL Test quaternion element equality
            %
            % ISEQUAL(Q1,Q2) is true if the Quaternions Q1 and Q2 are equal.
            %
            % Notes::
            % - Used by test suite verifyEqual() in addition to eq().
            % - Invokes eq() so respects double mapping for UnitQuaternion.
            %
            % See also Quaternion.eq.
            e = eq(q1, q2);
        end

        function e = eq(q1, q2)
            %EQ Test quaternion equality
            %
            % Q1 == Q2 is true if the Quaternions Q1 and Q2 are equal.
            %
            % Notes::
            % - Overloaded operator '=='.
            % - Equality means elementwise equality of Quaternion elements.
            % - If either, or both, of Q1 or Q2 are vectors, then the result is a vector.
            %   - if Q1 is a vector (1xN) then R is a vector (1xN) such that R(i) = Q1(i)==Q2.
            %   - if Q2 is a vector (1xN) then R is a vector (1xN) such that R(i) = Q1==Q2(i).
            %   - if both Q1 and Q2 are vectors (1xN) then R is a vector (1xN) such 
            %     that R(i) = Q1(i)==Q2(i).
            %
            % See also Quaternion.ne.
            if (numel(q1) == 1) && (numel(q2) == 1)
                e = sum(abs(q1.double - q2.double)) < 100*eps;
            elseif (numel(q1) >  1) && (numel(q2) == 1)
                e = zeros(1, numel(q1));
                for i=1:numel(q1)
                    e(i) = q1(i) == q2;
                end
            elseif (numel(q1) == 1) && (numel(q2) > 1)
                e = zeros(1, numel(q2));
                for i=1:numel(q2)
                    e(i) = q2(i) == q1;
                end
            elseif numel(q1) == numel(q2)
                e = zeros(1, numel(q1));
                for i=1:numel(q1)
                    e(i) = q1(i) == q2(i);
                end
            else
                error('SMTB:Quaternion:badarg', 'vectors not of same length');
            end
        end
        
        function e = ne(q1, q2)
            %NE Test quaternion inequality
            %
            % Q1 ~= Q2 is true if the Quaternions Q1 and Q2 are not equal.
            %
            % Notes::
            % - Overloaded operator '~='.
            % - If either, or both, of Q1 or Q2 are vectors, then the result is a vector.
            %   - if Q1 is a vector (1xN) then R is a vector (1xN) such that R(i) = Q1(i)~=Q2.
            %   - if Q2 is a vector (1xN) then R is a vector (1xN) such that R(i) = Q1~=Q2(i).
            %   - if both Q1 and Q2 are vectors (1xN) then R is a vector (1xN) such 
            %     that R(i) = Q1(i)~=Q2(i).
            %
            % See also Quaternion.eq.
            if (numel(q1) == 1) && (numel(q2) == 1)
                e = all( ne(q1.double, q2.double) );
            elseif (numel(q1) >  1) && (numel(q2) == 1)
                e = zeros(1, numel(q1));
                for i=1:numel(q1)
                    e(i) = q1(i) ~= q2;
                end
            elseif (numel(q1) == 1) && (numel(q2) > 1)
                e = zeros(1, numel(q2));
                for i=1:numel(q2)
                    e(i) = q2(i) ~= q1;
                end
            elseif numel(q1) == numel(q2)
                e = zeros(1, numel(q1));
                for i=1:numel(q1)
                    e(i) = q1(i) ~= q2(i);
                end
            else
                error('SMTB:Quaternion:badarg','vectors not of same length');
            end
        end
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TYPE CONVERSION METHODS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function s = char(q)
            %Quaternion.char Convert to string
            %
            % S = Q.char() is a compact string representation of the Quaternion's value
            % as a 4-tuple.  If Q is a vector then S has one line per element.
            %
            % Notes::
            % - The vector part is delimited by double angle brackets, to differentiate
            %   from a UnitQuaternion which is delimited by single angle brackets.
            %
            % See also UnitQuaternion.char.
            
            if length(q) > 1
                s = '';
                for qq = q;
                    s = char(s, char(qq));
                end
                return
            end
            
            function s = render(x)
                if isnumeric(x)
                    s = num2str(x);
                elseif isa(x, 'sym')
                    s = char(x);
                end
            end
                
            s = [render(q.s), ' << ' ...
                render(q.v(1)) ', ' render(q.v(2)) ', '   render(q.v(3)) ' >>'];
        end
                
        function v = double(q)
            %Quaternion.double Convert a quaternion to a 4-element vector
            %
            % V = Q.double() is a row vector (1x4) comprising the Quaternion elements,
            % scalar then vector, ie. V = [s vx vy vz].  If Q is a vector (1xN) of 
            % Quaternion objects then V is a matrix (Nx4) with rows corresponding to 
            % the quaternion elements.
            %
            
            for i=1:length(q)
                v(i,:) = [q(i).s q(i).v];
            end
        end
 
    end % methods
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STATIC FACTORY METHODS, ALTERNATIVE CONSTRUCTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static)
        
        function uq = new(varargin)
            %Quaternion.new Construct a new quaternion
            %
            % QN = Q.new() constructs a new Quaternion object.
            %
            % QN = Q.new([S, V1, V2, V3]) as above but specified directly by its 4 elements.
            %
            % QN = Q.new(S, V) as above but specified directly by the scalar S and vector
            % part V (1x3)
            %
            % Notes::
            % - Polymorphic with UnitQuaternion and RTBPose derived classes.
            uq = Quaternion(varargin{:});
        end
                
        function q = pure(v)
            %Quaternion.pure Construct a pure quaternion
            %
            % Q = Quaternion.pure(V) is a pure Quaternion formed from the vector V (1x3) and has
            % a zero scalar part.
            %

            
            if ~isvec(v)
                error('SMTB:Quaternion:badarg', 'must be a 3-vector');
            end
            q = Quaternion(0, v(:));
        end
    end  % static methods
end
