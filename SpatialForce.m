%SpatialForce Spatial force class
%
% Concrete subclass of SpatiallVec6 and SpatialF6 and represents the
% translational and rotational forces and torques acting on a rigid-body in 3D space.
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
%  SpatialForce           ^constructor invoked by subclasses
%  new                    construct new concrete class of same type
%  double                 ^convert to a 6xN double
%  char                   ^convert to string
%  display                ^display in human readable form
%
% Operators::
%  +          ^add spatial vectors of the same type
%  -          ^subtract spatial vectors of the same type
%  -          ^unary minus of spatial vectors
%  *          ^^^premultiplication by SE3 yields transformed SpatialForce
%  *          ^^^^premultiplication by Twist yields transformed SpatialForce
%
% Notes:
% - The implementation of methods indicated with ^ is inherited from SpatialVec6.
% - The implementation of methods indicated with ^^ is inherited from SpatialM6.
% - The implementation of methods indicated with ^^^ are implemented in RTBPose.
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
% See also SpatialVec6, SpatialF6, SpatialMomentum.

classdef SpatialForce < SpatialF6
    
    
    methods
        function n = new(a, val)
            n = SpatialForce(val);
        end
    end
end

