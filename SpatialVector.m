classdef SpatialVector < handle
    properties
        vw
    end
    
    methods
        function obj = SpatialVector(x)
            assert(isvec(x, 6), 'Must be a 6-vector');
            
            obj.vw = x(:);
        end
        
        function display(obj)
            disp( char(obj) );
        end

        function s = char(obj)
            s =sprintf('%s: [ %g %g %g | %g %g %g ]', class(obj), obj.vw);
        end
            
        function y = uminus(obj)
            y = obj.new(-obj.vw);
        end
        
        function y = plus(a, b)
            assert(strcmp(class(a), class(b)), 'can only add spatial vectors of same class')
            y = a.new(a.vw + b.vw);
        end
        
        function y = minus(a, b)
            assert(strcmp(class(a), class(b)), 'can only subtract spatial vectors of same class')
            v = a.new(a.vw - b.vw);
        end
                    
        function v = double(obj)
            v = obj.vw;
        end
    end
end
