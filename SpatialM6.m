%SpatialM6 Abstract spatial motion class
%
% Abstract superclass that represents spatial motion.  This class has two
% concrete subclasses:
%
%          SpatialVec6 (abstract handle class)
%            |
%            +--- SpatialM6 (abstract)
%            |     |
%            |     +---SpatialVelocity
%            |     +---SpatialAcceleration
%            |
%            +---SpatialF6 (abstract)
%                 |
%                 +---SpatialForce
%                 +---SpatialMomentum
%
% Methods::
%  SpatialM6     ^constructor invoked by subclasses
%  char          ^convert to string
%  cross         cross product
%  display       ^display in human readable form
%  double        ^convert to a 6xN double
%
% Operators::
%  +          ^add spatial vectors of the same type
%  -          ^subtract spatial vectors of the same type
%  -          ^unary minus of spatial vectors
%
% Notes:
%  - ^ is inherited from SpatialVec6.
%  - Subclass of the MATLAB handle class which means that pass by reference semantics
%    apply.
%  - Spatial vectors can be placed into arrays and indexed.
%
% References::
%
%  - Robot Dynamics Algorithms, R. Featherstone, volume 22,
%    Springer International Series in Engineering and Computer Science,
%    Springer, 1987.
%  - A beginner's guide to 6-d vectors (part 1), R. Featherstone, 
%    IEEE Robotics Automation Magazine, 17(3):83-94, Sep. 2010.
%
% See also SpatialForce, SpatialMomentum, SpatialInertia, SpatialM6.

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

classdef (Abstract) SpatialM6 < SpatialVec6
    methods
        

        function out = cross(obj, other)
            %SpatialM6.cross Spatial velocity cross product
            %
            % cross(V1, V2) is a SpatialAcceleration object where V1 and V2 are SpatialM6
            % subclass instances.
            %
            % cross(V, F) is a SpatialForce object where V1 is a SpatialM6
            % subclass instances and F is a SpatialForce subclass instance.
            %
            % Notes::
            %  - The first form is Featherstone's "x" operator.
            %  - The second form is Featherstone's "x*" operator.

            v = obj.vw;
            % vcross = [ skew(w) skew(v); zeros(3,3) skew(w) ]
            
            
            vcross = [
                0    -v(6)  v(5)     0    -v(3)  v(2)
                v(6)  0    -v(4)     v(3)  0    -v(1)
                -v(5)  v(4)  0       -v(2)  v(1)  0
                0     0     0        0    -v(6)  v(5)
                0     0     0       v(6)  0    -v(4)
                0     0     0      -v(5)  v(4)  0
                ];
            if isa(other, 'SpatialVelocity')
                out = SpatialAcceleration( vcross * other.vw );
            elseif isa(other, 'SpatialF6')
                vcross = -vcross';      % x* operator
                out = SpatialForce( vcross * other.vw );
            else
                error('type mismatch')
            end
            
        end
        
    end
end

