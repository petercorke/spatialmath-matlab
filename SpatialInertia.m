%SpatialInertia Spatial inertia class
%
% Concrete class representing spatial inertia.
%
% Methods::
%  SpatialInertia   constructor 
%  char             convert to string
%  display          display in human readable form
%  double           convert to a 6xN double
%
% Operators::
%  +          plus: add spatial inertia of connected bodies
%  *          mtimes: compute force or momentum 
%
% Notes::
%  - Subclass of the MATLAB handle class which means that pass by reference semantics
%    apply.
%  - Spatial inertias can be placed into arrays and indexed.
%
% References::
%
%  - Robot Dynamics Algorithms, R. Featherstone, volume 22,
%    Springer International Series in Engineering and Computer Science,
%    Springer, 1987.
%  - A beginner's guide to 6-d vectors (part 1), R. Featherstone, 
%    IEEE Robotics Automation Magazine, 17(3):83-94, Sep. 2010.
%
% See also SpatialM6, SpatialF6, SpatialVelocity, SpatialAcceleration, SpatialForce,
% SpatialMomentum.

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

classdef SpatialInertia < handle
    
    properties
        I
    end
    
    methods
        function obj = SpatialInertia(m, c, I)
            %SpatialInertia.SpatialInertia Constructor
            %
            % SI = SpatialInertia(M, C, I) is a spatial inertia object for a rigid-body
            % with mass M, centre of mass at C relative to the link frame, and an
            % inertia matrix (3x3) about the centre of mass.
            %
            % SI = SpatialInertia(I) is a spatial inertia object with a value equal
            % to I (6x6).
            switch nargin
                case 0
                    obj.I = zeros(3,3);
                case 1
                    assert(all(size(m) == [6 6]), 'Must pass a 6x6 matrix');
                    
                    obj.I = m;
                case 3
                    C = skew(c);
                    
                    obj.I = [
                        m*eye(3) m*C'
                        m*C I+m*C*C'
                        ];
            end
        end
        
        function display(obj)
            %SpatialInertia.display Display parameters
            %
            % SI.display() displays the spatial inertia parameters in compact format.
            % If SI is an array of spatial inertia objects it displays them in a vertical
            % list.
            %
            % Notes::
            %  - This method is invoked implicitly at the command line when the result
            %    of an expression is a spatial inerita object and the command has
            %    no trailing semicolon.
            %
            % See also SpatialInertia.char.
            
            loose = strcmp( get(0, 'FormatSpacing'), 'loose');
            if loose
                disp(' ');
            end
            disp([inputname(1), ' = '])
            disp( char(obj) )
        end

        function s = char(obj, flag)
             %SpatialInertia.char Convert to string
             %
             % s = SI.char() is a string showing spatial inertia parameters in a 
             % compact format.
             % If SI is an array of spatial inertia objects return a string with the
             % inertia values in a vertical list.
             %
             % See also SpatialInertia.display.

             if numel(obj) == 1
                 if nargin == 1 || flag == 1
                    s = sprintf('%s:', class(obj));
                 else
                     s = '';
                 end
                 m = num2str(obj.I);
                 for line = m'
                     s = char(s, ['    ' line']);
                 end
             else
                 s = char( obj(1) );
                 
                 for i = 2:numel(obj)
                     s = char(s, ' ');
                     s = char(s, char(obj(i), 0) );
                 end
             end
        end
        
        function v = plus(a,b)
            %SpatialInertia.plus Addition operator
            %
            % SI1 + SI2 is the SpatialInertia of a composite body when bodies with SpatialInertia SI1 and SI2
            % are connected.
            %
            assert(isa(b, 'SpatialInertia'), 'spatial inertia can only be added to spatial inertia')
            v = SpatialInertia( a.I + b.I );
        end
        
        function v = mtimes(a,b)
            %SpatialInertia.times Multiplication operator
            %
            % SI * A is the SpatialForce required for a body with SpatialInertia SI to accelerate with 
            % the SpatialAcceleration A.
            %
            % SI * V is the SpatialMomemtum of a body with SpatialInertia SI and SpatialVelocity V.
            %
            % Notes::
            % - These products must be written in this order, A*SI and V*SI are not defined.

            if isa(b, 'SpatialAcceleration')
                    v = SpatialForce(a.I * b.vw);   % F = ma
            elseif isa(b, 'SpatialVelocity')
                    % crf(v(i).vw)*model.I(i).I*v(i).vw;
                    %v = Wrench( a.cross() * I.I * a.vw );
                    v = SpatialMomentum(a.I * b.vw);
            else
                error( 'bad postmultiply operands for Inertia *'); % M = mv
            end
            
        end
        
        function v = double(obj)
            %SpatialInertia.double Convert to matrix
            %
            % double(V) is a native matrix (6x6) with the value of the spatial inertia.
            % If V is an array (1xN) the result is a matrix (6x6xN).
            v = reshape( [obj.I], 6, 6, []);
        end
    end
end