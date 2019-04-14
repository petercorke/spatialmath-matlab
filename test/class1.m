classdef class1 < handle

    properties
        x
    end


    methods
        function obj = class1(x)
            obj.x = x;
        end

        function v = thing(x)
            97
            fprintf('  thing -> %d\n', length(x))
            v = 12;
            98
        end

        function disp(obj)
            87
            fprintf('value is %f\n', obj.x);
            88
        end
    end
end
