classdef SpatialInertia < handle
    
    properties
        I
    end
    
    methods
        function obj = SpatialInertia(m, c, I)
            if nargin < 3
                I = zeros(3,3);
            end
            C = skew(c);
            
            obj.I = [
                m*eye(3) m*C'
                m*C I+m*C*C'
                ];
        end
        
        function display(obj)
            disp( char(obj) );
        end

        function s = char(obj)
            s =sprintf('%s:', class(obj));
            m = num2str(obj.I);
            for line = m'
                s = strvcat(s, line');
            end
        end
        
        function v = mtimes(a,b)
            % Inertia * a -> Wrench
            % Inertia * v -> Wrench

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
    end
end