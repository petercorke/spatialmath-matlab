%SpatialVelocity Spatial velocity class
%
% Concrete subclass of SpatiallVec6 and SpatialM6 and represents the
% translational and rotational velocity of a rigid-body moving in 3D space.
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
%  SpatialVelocity        ^constructor invoked by subclasses
%  new                    construct new concrete class of same type
%  double                 ^convert to a 6xN double
%  cross                  ^^cross product
%  char                   ^convert to string
%  display                ^display in human readable form
%
% Operators::
%
%  +          ^add spatial vectors of the same type
%  -          ^subtract spatial vectors of the same type
%  -          ^unary minus of spatial vectors
%  *          ^^^premultiplication by SpatialInertia yields SpatialMomentum
%  *          ^^^^premultiplication by Twist yields transformed SpatialVelocity
%
% Notes:
% - The implementation of methods indicated with ^ is inherited from SpatialVec6.
% - The implementation of methods indicated with ^^ is inherited from SpatialM6.
% - The implementation of methods indicated with ^^^ are implemented in SpatialInertia.
% - The implementation of methods indicated with ^^^^ are implemented in Twist.
%
% References::
%
% - Robot Dynamics Algorithms, R. Featherstone, volume 22,
%   Springer International Series in Engineering and Computer Science,
%   Springer, 1987.
% - A beginner?s guide to 6-d vectors (part 1), R. Featherstone, 
%   IEEE Robotics Automation Magazine, 17(3):83?94, Sep. 2010.
%
% See also SpatialVec6, SpatialM6, SpatialAcceleration, SpatialInertia, SpatialMomentum.

classdef SpatialVelocity < SpatialM6
    methods

        function n = new(a, val)
            n = SpatialVelocity(val);
        end
    end
end
