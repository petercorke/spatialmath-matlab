%SpatialM6 Abstract spatial motion class
%
% Abstract superclass that represents spatial motion.  This class has two
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
%  SpatialM6     ^constructor invoked by subclasses
%  double        ^convert to a 6xN double
%  char          ^convert to string
%  display       ^display in human readable form
%
% Operators::
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

classdef (Abstract) SpatialM6 < SpatialVec6
    methods
        
        %SpatialM6.cross Spatial velocity cross product
        %
        % cross(V, V) is a SpatialAcceleration object.
        %
        % cross(V, F) is a SpatialForce object where F is a subclass of
        % SpatialForce.
        %
        % Notes::
        % - The first form is Featherstone's "x" operator.
        % - The second form is Featherstone's "x*" operator.
        function out = cross(obj, other)
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

