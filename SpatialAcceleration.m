%SpatialAcceleration Spatial acceleration class
%
% Concrete subclass of SpatiallVec6 and SpatialM6 and represents the
% translational and rotational acceleration of a rigid-body moving in 3D space.
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
%  SpatialAcceleration    ^constructor invoked by subclasses
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
%  *          ^^^premultiplication by SpatialInertia yields SpatialForce
%  *          ^^^^premultiplication by Twist yields transformed SpatialAcceleration
%
%
% Notes:
% - The implementation of methods indicated with ^ is inherited from SpatialVec6.
% - The implementation of methods indicated with ^^ is inherited from SpatialM6.
% - The implementation of methods indicated with ^^^ are implemented in SpatialInertia.
% - The implementation of methods indicated with ^^^^ are implemented in Twist.

classdef SpatialAcceleration < SpatialM6
    
    methods
        function n = new(a, val)
            n = SpatialAcceleration(val);
        end
    end
end
