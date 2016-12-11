

% plot_ellipse
%     2d outline, filled case
%     3d outlien, filled case
%     with LS or edgecolor, color options etc.

function tests = TransformationsTest
  tests = functiontests(localfunctions);
end

% plot_ellipse
function ellipse2d_test(testCase)
    clf
    
    E = diag([9 4]);
    
    plot_ellipse(E)
    plot_ellipse(E, [4 0], 'r--');
    plot_ellipse(E, [0 4], 'edgecolor', 'g');
    plot_ellipse(E, [4 4], 'fillcolor', 'g');
    plot_ellipse(E, [0 8 0.5], 'edgecolor', 'r', 'fillcolor', 'c', 'alpha', 0.5, 'LineWidth', 3);
    plot_ellipse(E, [4 8], 'b--', 'LineWidth', 3);

    axis equal
end

function ellipse2d_animate_test(testCase)
    clf
    axis([-4 4 -4 4]);
    
    E = diag([9 4]);
    
    h = plot_ellipse(E, 'g');
    for x = circle([0 0], 1)
        plot_ellipse(E, x, 'alter', h)
        pause(0.1)
    end
 
        clf
    axis([-4 4 -4 4]);
    h = plot_ellipse(E, [0 0 0.5], 'edgecolor', 'r', 'fillcolor', 'c', 'LineWidth', 3);

    for x = circle([0 0], 1)
        plot_ellipse(E, x, 'alter', h)
        pause(0.1)
    end

end

% plot_ellipse
function ellipse3d_test(testCase)
    clf
    
    E = diag([9 4 6]);
    
    plot_ellipse(E)
    pause
    
    clf
    plot_ellipse(E, 'edgecolor', 'g');
    pause
    
    clf
    plot_ellipse(E, 'fillcolor', 'g');
    pause
    
    clf
    plot_ellipse(E, 'fillcolor', 'g', 'shadow');
    pause
    
    clf
    plot_ellipse(E, 'fillcolor', 'g', 'edgecolor', 'r', 'LineWidth', 2);
    pause
    
    plot_ellipse(E, [0 8], 'edgecolor', 'r', 'fillcolor', 'c');
    plot_ellipse(E, [4 8], 'LineWidth', 3, 'MarkerStyle', '+');

    axis equal
end

function ellipse3d_animate_test(testCase)
    clf
    axis([-4 4 -4 4 -4 4]);
    
    E = diag([9 4 6]);
    
    h = plot_ellipse(E, 'g');
    for x = circle([0 0], 1)
        plot_ellipse(E, [x; 0], 'alter', h)
        pause(0.1)
    end
end