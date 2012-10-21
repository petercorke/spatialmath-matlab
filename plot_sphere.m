%PLOT_SPHERE Plot spheres
%
% PLOT_SPHERE(C, R, COLOR) add spheres to the current figure.  C is the 
% centre of the sphere and if its a 3xN matrix then N spheres are drawn 
% with centres as per the columns.  R is the radius and COLOR is a Matlab 
% color spec, either a letter or 3-vector.
%
% H = PLOT_SPHERE(C, R, COLOR) as above but returns the handle(s) for the
% spheres.
%
% H = PLOT_SPHERE(C, R, COLOR, ALPHA) as above but ALPHA specifies the opacity
% of the sphere were 0 is transparant and 1 is opaque.  The default is 1.
%
% Example::
% Create four spheres
%         plot_sphere( mkgrid(2, 1), .2, 'b')
% and now turn on a full lighting model
%         lighting gouraud
%         light
%
% NOTES::
% - The sphere is always added, irrespective of figure hold state.
% - The number of vertices to draw the sphere is hardwired.

% TODO
% inconsistant call format compared to other plot_xxx functions.

function h = plot_sphere(c, r, varargin)

    opt.color = 'b';
    opt.alpha = 1;
    opt.mesh = 'none';

    [opt,args] = tb_optparse(opt, varargin);
    
    % backward compatibility with RVC
    if length(args) > 0
        opt.color = args{1};
    end
    if length(args) > 1
        opt.alpha = args{2};
    end
    
    daspect([1 1 1])
    hold_on = ishold
    hold on
    [xs,ys,zs] = sphere(40);

    if isvec(c,3)
        c = c(:);
    end
    if size(r) == 1
        r = r * ones(numcols(c),1);
    end

    if nargin < 4
        alpha = 1
    end

    % transform the sphere
    for i=1:numcols(c)
        x = r(i)*xs + c(1,i);
        y = r(i)*ys + c(2,i);
        z = r(i)*zs + c(3,i);
                
        % the following displays a nice smooth sphere with glint!
        h = surf(x,y,z, 'FaceColor', opt.color, 'EdgeColor', opt.mesh, 'FaceAlpha', opt.alpha);
        % camera patches disappear when shading interp is on
        %h = surfl(x,y,z)
    end
    %lighting gouraud
    %light
    if ~hold_on
        hold off
    end
