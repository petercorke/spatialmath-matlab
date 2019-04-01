%SpatialVelocity Spatial velocity class
%
% Concrete subclass of SpatialM6 and represents the
% translational and rotational velocity of a rigid-body moving in 3D space.
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
%  SpatialVelocity    ^constructor invoked by subclasses
%  char               ^convert to string
%  cross              ^^cross product
%  display            ^display in human readable form
%  double             ^convert to a 6xN double
%  new                construct new concrete class of same type
%
% Operators::
%
%  +      ^add spatial vectors of the same type
%  -      ^subtract spatial vectors of the same type
%  -      ^unary minus of spatial vectors
%  *      ^^^premultiplication by SpatialInertia yields SpatialMomentum
%  *      ^^^^premultiplication by Twist yields transformed SpatialVelocity
%
% Notes:
% - ^ is inherited from SpatialVec6.
% - ^^ is inherited from SpatialM6.
% - ^^^ are implemented in SpatialInertia.
% - ^^^^ are implemented in Twist.
%
% References::
%
%  - Robot Dynamics Algorithms, R. Featherstone, volume 22,
%    Springer International Series in Engineering and Computer Science,
%    Springer, 1987.
%  - A beginner's guide to 6-d vectors (part 1), R. Featherstone, 
%    IEEE Robotics Automation Magazine, 17(3):83-94, Sep. 2010.
%
% See also SpatialVec6, SpatialM6, SpatialAcceleration, SpatialInertia, SpatialMomentum.

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

classdef SpatialVelocity < SpatialM6
    methods

        function n = new(a, val)
            %SpatialVelocity.new Construct a new object of the same type
            %
            % A2 = A.new(X) creates a new object of the same type as A, with the value 
            % X (6x1).
            %
            % Notes::
            %  - Serves as a dynamic constructor.
            %  - This method is polymorphic across all SpatialVec6 derived classes, and
            %    allows easy creation of a new object of the same class as an existing
            %    one without needing to explicitly determine its type.
            n = SpatialVelocity(val);
        end
    end
end
