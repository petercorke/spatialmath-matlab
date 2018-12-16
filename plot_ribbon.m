clf
trplot( SE3 )

plot_ribbon2( SE3 )


%   ^
%   |                            | \
%   |  ^  +----------------------|  \
%   |  |  |                          .
%   |  v  +----------------------|  /
%   v  w1                        | /
%   w2    <---------- d --------->
%         <--------- wrap ----------->
function plot_ribbon2(pose)
    radius = 0.5;
    d = 0.90;
    w1 = 0.2;
    w2 = 0.4;
    N = 100;
    Nd = floor(N*d);
    wrap = 0.8;
    phase = -.3;
    
    % compute canonic arrow about the z-axis and centered at origin
    
    
    theta = linspace(0, wrap*2*pi, N) + phase + wrap*2*pi/2;
    % replicate the Nd'th point, this means that
    % x(Nd) = x(Nd+1)
    % y(Nd) = y(Nd+1)
    % which makes the base of the arrowhead normal to the shaft
    theta = [theta(1:Nd) theta(Nd) theta(Nd+1:end)];
    x = radius * cos(theta);
    y = radius * sin(theta);
    
    % compute the width of the arrow, varies along length
    T = [0 w1/2; Nd-1 w1/2; Nd w2/2; N-1 0];
    z = interp1(T(:,1), T(:,2), [0:N], 'linear');
    
    
    % build the mesh matrices
    X = [x; x];
    Y = [y; y];
    Z = [-z; z];
    
    z1 = ones(size(Z));
    C = cat(3, z1*1, z1*0, z1*0);
    hg = hgtransform
    hold on
    surf(X, Y, Z, C, 'Parent', hg, 'EdgeColor', 'None', 'FaceColor', [1 0 0])
    hold off
    
    hg.Matrix = pose.T;
end