%SpatialV6 Abstract spatial 6-vector class
%
% Abstract superclass for spatial vector functionality.  This class has two
% abstract subclasses, which each have concrete subclasses:
%
%     SpatialVec6 (abstract handle class)
%        |
%        +--- SpatialM6 (abstract)
%        |     |
%        |     +---SpatialVelocity
%        |     +---SpatialAcceleration
%        |
%        +---SpatialF6 (abstract)
%             |
%             +---SpatialForce
%             +---SpatialMomentum
%
% Methods::
%  SpatialV6     constructor invoked by subclasses
%  double        convert to a 6xN double
%  char          convert to string
%  display       display in human readable form
%
% Operators::
%
%  +          add spatial vectors of the same type
%  -          subtract spatial vectors of the same type
%  -          unary minus of spatial vectors
%
% Notes::
% - Subclass of the MATLAB handle class which means that pass by reference semantics
%   apply.
% - Spatial vectors can be placed into arrays and indexed.
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
% SpatialMomentum, SpatialInertia.
%

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
classdef (Abstract) SpatialVec6 < handle
    properties
        vw  % the 6-vector: translation elements, then rotation
    end
    
    methods
        function obj = SpatialVec6(x)
            %SpatialVec6.SpatialVec6 Constructor
            %
            % SpatiaVecXXX(V) is a spatial vector of type SpatiaVecXXX with a value
            % from V (6x1).  If V (6xN) then an (Nx1) array of spatial vectors is
            % returned.
            %
            % This constructor is inherited by all the concrete subclasses.
            %
            % See also SpatialVelocity, SpatialAcceleration, SpatialForce, SpatialMomentum.
            
            assert((numel(x)==6) || (size(x,1) == 6), 'Must be a 6-vector or matrix with 6 rows');
            if length(x) == 6
                obj.vw = x(:);
            else
                n = size(x, 2);
                for i=1:n
                    obj(i) = obj.new(x(:,i));
                end
            end
            
        end
        
        function display(obj)
            %SpatialVec6.display Display parameters
            %
            % V.display() displays the spatial vector parameters in compact single line format.
            % If V is an array of spatial vector objects it displays one per line.
            %
            % Notes::
            %  - This method is invoked implicitly at the command line when the result
            %    of an expression is a serial vector subclass object and the command has
            %    no trailing semicolon.
            %
            % See also SpatialVec6.char.
            
            loose = strcmp( get(0, 'FormatSpacing'), 'loose');
            if loose
                disp(' ');
            end
            disp([inputname(1), ' = '])
            disp( char(obj) )
        end

        function s = char(obj)
             %SpatialVec6.char Convert to string
             %
             % s = V.char() is a string showing spatial vector parameters in a 
             % compact single line format.
             % If V is an array of spatial vector objects return a string with one
             % line per element.
             %
             % See also SpatialVec6.display.
            if numel(obj) == 1
                s = sprintf('%s: [ %g %g %g | %g %g %g ]', class(obj), obj.vw);
            else
                s = char( obj(1) );
                indent = repmat(' ', 1, length(class(obj)));
                for i = 2:numel(obj)
                    s = char(s, sprintf('%s  [ %g %g %g | %g %g %g ]', indent, obj(i).vw));
                end
            end
        end
            
        function y = uminus(obj)
            %SpatialVec6.uminus Unary minus operator
            %
            % -V is a spatial vector of the same type as V whose value is 
            % the negative of V.  If V is an array V (1xN) then the result 
            % is an array (1xN).
            %
            % See also SpatialVec6.minus, SpatialVec6.plus.
            for i=1:numel(obj)
                y(i) = obj.new(-obj(i).vw);
            end
        end
        
        function y = plus(a, b)
            %SpatialVec6.plus Addition operator
            %
            % V1 + V2 is a spatial vector of the same type as V1 and V2 whose value is 
            % the sum of V1 and V2.  If both are arrays of spatial vectors V1 (1xN) and
            % V2 (1xN) the result is an array (1xN).
            %
            % See also SpatialVec6.minus.
            
            assert(strcmp(class(a), class(b)), 'can only add spatial vectors of same type')
            assert(numel(a) == numel(b), 'can only add equal length arrays of spatial vectors');
            for i=1:numel(a)
                y(i) = a.new(a(i).vw + b(i).vw);
            end
        end
        
        function y = minus(a, b)
            %SpatialVec6.minus Subtraction operator
            %
            % V1 - V2 is a spatial vector of the same type as V1 and V2 whose value is 
            % the difference of V1 and V2.  If both are arrays of spatial vectors V1 (1xN) and
            % V2 (1xN) the result is an array (1xN).
            %
            % See also SpatialVec6.uminus, SpatialVec6.plus.
            assert(strcmp(class(a), class(b)), 'can only subtract spatial vectors of same class')
            assert(numel(a) == numel(b), 'can only subtract equal length arrays of spatial vectors');
            for i=1:numel(a)
                y(i) = a.new(a(i).vw - b(i).vw);
            end
        end
                    
        function v = double(obj)
            %SpatialVec6.dpuble Convert to matrix
            %
            % double(V) is a native matrix (6x1) with the value of the spatial vector.
            % If V is an array (1xN) the result is a matrix (6xN).
            v = [obj.vw];
        end
    end
end
