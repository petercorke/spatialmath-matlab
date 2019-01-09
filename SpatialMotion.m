classdef SpatialMotion < SpatialVector
    methods
        
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
            elseif isa(other, 'SpatialForce')
                vcross = -vcross';      % x* operator
                out = SpatialForce( vcross * other.vw );
            else
                error('type mismatch')
            end
            
            
            
        end
        
    end
end

