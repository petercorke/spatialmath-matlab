classdef Animate < handle
    properties
        frame
        dir
    end
    
    methods
        function a = Animate(name)
            a.frame = 0;
            a.dir = name;
            mkdir(name);
            delete( fullfile(name, '*.png') );
            
        end
        
        function add(a, fh)
            if nargin < 2
                fh = gcf;
            end
            
            print('-dpng', fullfile(a.dir, sprintf('%04d.png', a.frame)));
            a.frame = a.frame + 1;
        end
    end
end