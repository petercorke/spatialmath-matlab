%ANIMATE Create an animation
%
% Helper class for creating animations.  Creates a movie file or saves snapshots 
% of a figure as individual PNG format frames numbered 0000.png, 0001.png and so
% on.
%
% Example::
%
%         anim = Animate('movie.mp4');
%          for i=1:100
%              plot(...);
%              anim.add();
%          end
%          anim.close();
%
% will save the frames in an MP4 movie file using VideoWriter.
%
% Alternatively::
%
%          anim = Animate('movie');
%          for i=1:100
%              plot(...);
%              anim.add();
%          end
%          anim.close();
%
% To convert the image files to a movie you could use a tool like ffmpeg
%           ffmpeg -r 10 -i movie/%04d.png out.mp4
%


% Copyright (C) 1993-2017, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for MATLAB (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.
%
% http://www.petercorke.com


classdef Animate < handle
    properties
        frame
        dir
        resolution
        video
    end
    
    methods
        function a = Animate(name, varargin)
            %ANIMATE.ANIMATE Create an animation class
            %
            % A = ANIMATE(NAME, OPTIONS) initializes an animation, and creates
            % a movie file or a folder holding individual frames.
            %
            % A = ANIMATE({NAME, OPTIONS}) as above but arguments are passed as a cell array,
            % which allows a single argument to a higher-level option like 'movie',M to express
            % options as well as filename.
            %
            % Options::
            % 'resolution',R     Set the resolution of the saved image to R pixels per inch.
            % 'profile',P        Create an MP4 file directly, see VideoWriter
            % 'fps',F            Frame rate (default 30)
            %
            % Notes::
            % - if a profile is given a movie is created, see VideoWriter for allowable
            %   profiles.
            % - if the file has an extension a movie is created.
            % - otherwise a folder full of frames is created.
            %
            % See also VideoWriter.
        
            if isempty(name)
                % we're not animating
                a.dir = [];
            else
                opt.resolution = [];
                opt.profile = [];
                opt.fps = 30;

                if iscell(name) && length(varargin) == 0
                    varargin = name(2:end);
                    name = name{1};
                end
                [opt,args] = tb_optparse(opt, varargin);
                a.resolution = opt.resolution;
                a.frame = 0;

                [p,f,e] = fileparts(name);
                
                if ~isempty(opt.profile)
                    % create a video with this profile
                    a.video = VideoWriter(name, a.profile, args{:});
                    fprintf('saving video --> %s with profile ''%s''\n', name, a.profile);

                elseif ~isempty(e)
                    % an extension was given
                    switch (e)
                        case {'.mp4', '.m4v'},  profile = 'MPEG-4';
                        case '.mj2',  profile = 'Motion JPEG 2000';
                        case '.avi', profile = 'Motion JPEG AVI';
                    end
                    fprintf('Animate: saving video --> %s with profile ''%s''\n', name, profile);
                    a.video = VideoWriter(name, profile, args{:});
                    a.video.FrameRate = opt.fps;
                    a.video.Quality = 95;
                    open(a.video);
                else
                    % create a folder to hold the frames
                    a.dir = name;
                    mkdir(name);
                    
                    % clear out old frames
                    delete( fullfile(name, '*.png') );
                    fprintf('saving frames --> %s\n', name);

                end
            end
            
        end
        
        function add(a, fh)
            %ANIMATE.ADD Adds current plot to the animation
            %
            % A.ADD() adds the current figure in PNG format to the animation
            % folder with a unique sequential filename.
            %
            % A.ADD(FIG) as above but captures the figure FIG.
            %
            % See also print.

            if isempty(a.dir) && isempty(a.video)
                return;
            end
            
            if nargin < 2
                fh = gcf;
            end
            
            if isempty(a.video)
                if isempty(a.resolution)
                    print(fh, '-dpng', fullfile(a.dir, sprintf('%04d.png', a.frame)));
                else
                    print(fh, '-dpng', sprintf('-r%d', a.resolution), fullfile(a.dir, sprintf('%04d.png', a.frame)));
                end
            else
                im = frame2im( getframe(fh) );  % get the frame
                
                % crop so that height/width are multiples of two, by default MATLAB pads
                % with black which gives lines at the edge
                w = numcols(im); h = numrows(im);
                w = floor(w/2)*2; h = floor(h/2)*2;
                im = im(1:h,1:w,:);
                
                % add the frame to the movie
                writeVideo(a.video, im)
            end
            a.frame = a.frame + 1;
        end
        
        function out = close(a)
            %ANIMATE.CLOSE  Closes the animation
            %
            % A.CLOSE() closes the video file.
            %
            % 
            if ~isempty(a.video)
                if nargout > 0
                    out = char(a.video);
                end
                close(a.video);
            end
        end
    end
end
