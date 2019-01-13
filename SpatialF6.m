%SpatialF6 Abstract spatial force class
%
% Abstract superclass that represents spatial force.  This class has two
% concrete subclasses:
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
%  SpatialF6     ^constructor invoked by subclasses
%  double        ^convert to a 6xN double
%  char          ^convert to string
%  display       ^display in human readable form
%
% Operators::
%
%  +          ^add spatial vectors of the same type
%  -          ^subtract spatial vectors of the same type
%  -          ^unary minus of spatial vectors
%
% Notes:
% - The implementation of methods indicated with ^ is inherited from SpatialVec6.
% - Subclass of the MATLAB handle class which means that pass by reference semantics
%   apply.
% - Spatial vectors can be placed into arrays and indexed.
%
% References::
%
% - Robot Dynamics Algorithms, R. Featherstone, volume 22,
%   Springer International Series in Engineering and Computer Science,
%   Springer, 1987.
% - A beginner?s guide to 6-d vectors (part 1), R. Featherstone, 
%   IEEE Robotics Automation Magazine, 17(3):83?94, Sep. 2010.
%
% See also SpatialForce, SpatialMomentum, SpatialInertia, SpatialM6.

classdef (Abstract) SpatialF6 < SpatialVec6
end

