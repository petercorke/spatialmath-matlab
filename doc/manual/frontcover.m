function plot_arrow2(p1, p2, varargin)
    p1 = [0 0 0];
    p2 = [0 0 2];
    clf

    opt.color = 'b';
    opt.edgecol = 'none';
    opt.shaftdiam = 0.2;
    opt.headwidthscale = 2;
    opt.headlengthscale = 0.2;
    
        l = norm(p2-p1)
    
    h1 = arrow(3, opt, 'r'); h1.Matrix = troty(pi/2);
    h2 = arrow(3, opt, 'g'); h2.Matrix = trotx(-pi/2);
    h3 = arrow(3, opt, 'b');
    
    ht = hgtransform()
    h1.Parent = ht;
    h2.Parent = ht;
    h3.Parent = ht;
    
    Tf = transl(7.5,7, 8) * rpy2tr(0.6, 0.8, 1.2);
    N = 50
    T = trinterp(Tf, N);
    hold on
    for i=1:N-1
        trplot(T(:,:,i), 'notext', 'color', [1 1 1]*0.3, 'perspective');
    end
    ht.Matrix = T(:,:,end);
            view(30, -13)
    %shading interp
    lighting gouraud
    light
    axis equal
            axis([0 10 0 10 0 10])
            grid on
            xlabel('X', 'FontSize', 16);
            ylabel('Y', 'FontSize', 16);
            zlabel('Z', 'FontSize', 16);
            view(158, 21)

end

function hg = arrow(l, opt, col)
    
    if nargin > 2
        opt.color = col;
    end
    
    hg = hgtransform()
    
    hl = opt.headlengthscale*l
    hw = opt.shaftdiam*opt.headwidthscale
    
    % shaft
    [x,y,z] = cylinder(opt.shaftdiam/2);
    z = z*(l-hl);  % shaft length
    surf(x,y,z, ones(size(z)), 'FaceColor', opt.color, 'EdgeColor', opt.edgecol, 'parent', hg);
    
    % arrow head
    [x,y,z] = cylinder([hw 0]/2);
    z = z*hl + (l-hl);
    surf(x,y,z, ones(size(z)), 'FaceColor', opt.color, 'EdgeColor', opt.edgecol, 'parent', hg);
    
    patch('XData', x(1,:), 'YData', y(1,:), 'ZData', z(1,:), 'FaceColor', opt.color, 'EdgeColor', opt.edgecol, 'parent', hg)
end