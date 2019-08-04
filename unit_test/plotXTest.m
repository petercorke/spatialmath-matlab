function tests = plotXTest
  tests = functiontests(localfunctions);
  
  clc
end

%     2d outline, filled case
%     3d outlien, filled case
%     with LS or edgecolor, color options etc.

function teardownOnce(tc)
    close all
end

function plotpoint_test(tc)
    
    % simple
    points = rand(2,5);
    clf; plot_point(points);
    tc.verifyEqual(length(get(gca, 'Children')), 5);
    tc.verifyEqual(length(findobj(gca, 'Type', 'line')), 5);
    lines = findobj(gca, 'Type', 'line');
    for i=1:5
        tc.verifyEqual(length(lines(i).XData), 1);
        tc.verifyEqual(length(lines(i).YData), 1);
        tc.verifyEqual(lines(i).LineStyle, 'none');
        tc.verifyEqual(lines(i).MarkerSize, 6);
        tc.verifyEqual(lines(i).MarkerFaceColor, 'none');
        tc.verifyEqual(lines(i).Marker, 'square');
    end
    
    % markers specified
    clf; plot_point(points, 'rd');
    tc.verifyEqual(length(get(gca, 'Children')), 5);
    tc.verifyEqual(length(findobj(gca, 'Type', 'line')), 5);
    lines = findobj(gca, 'Type', 'line');
    for i=1:5
        tc.verifyEqual(lines(i).LineStyle, 'none');
        tc.verifyEqual(lines(i).MarkerSize, 6);
        tc.verifyEqual(lines(i).MarkerFaceColor, 'none');
        tc.verifyEqual(lines(i).Marker, 'diamond');
        tc.verifyEqual(lines(i).Color, [1 0 0]);
    end
    
    % markers specified solid
    clf; plot_point(points, 'rd', 'solid');
    tc.verifyEqual(length(get(gca, 'Children')), 5);
    tc.verifyEqual(length(findobj(gca, 'Type', 'line')), 5);
    lines = findobj(gca, 'Type', 'line');
    for i=1:5
        tc.verifyEqual(lines(i).LineStyle, 'none');
        tc.verifyEqual(lines(i).MarkerSize, 6);
        tc.verifyEqual(lines(i).MarkerFaceColor, [1 0 0]);
        tc.verifyEqual(lines(i).Marker, 'diamond');
        tc.verifyEqual(lines(i).Color, [1 0 0]);
    end

    % sequential labels
    clf; plot_point(points, 'sequence');
    tc.verifyEqual(length(get(gca, 'Children')), 10);
    tc.verifyEqual(length(findobj(gca, 'Type', 'text')), 5);
    labels = findobj(gca, 'Type', 'text');
    set = [1:5];
    for i=1:5
        set = setdiff(set, str2num(labels(i).String));
    end
    tc.verifyEmpty(set, 'Not all labels found');
    
    % specified labels
    L = {'A', 'B', 'C', 'D', 'E'};
    clf; plot_point(points, 'label', L);
    tc.verifyEqual(length(get(gca, 'Children')), 10);
    tc.verifyEqual(length(findobj(gca, 'Type', 'text')), 5);
    labels = findobj(gca, 'Type', 'text');
    set = [1:5];
    for i=1:5
        set = setdiff(set, double(strip(labels(i).String))-'A'+1);
    end
    tc.verifyEmpty(set, 'Not all labels found');
    
    % printf labels
    clf; plot_point(points, 'printf', {'label=%d', 21:25});
    tc.verifyEqual(length(get(gca, 'Children')), 10);
    tc.verifyEqual(length(findobj(gca, 'Type', 'text')), 5);
    labels = findobj(gca, 'Type', 'text');
    set = [21:25];
    for i=1:5
        l = strip(labels(i).String);
        tc.verifyEqual(l(1:6), 'label=');
        set = setdiff(set, str2num(l(7:end)));
    end
    tc.verifyEmpty(set, 'Not all labels found');
    
    % specify font size
    clf; plot_point(points, 'sequence', 'textsize', 23);
    tc.verifyEqual(length(get(gca, 'Children')), 10);
    
    tc.verifyEqual(length(findobj(gca, 'Type', 'line')), 5);
    tc.verifyEqual(length(findobj(gca, 'Type', 'text')), 5);
    labels = findobj(gca, 'Type', 'text');
    for i=1:5
        tc.verifyEqual(labels(i).FontSize, 23);
        tc.verifyEqual(labels(i).FontWeight, 'normal');
    end
    
    % specify font color
    clf; plot_point(points, 'sequence', 'textcolor', 'g');
    tc.verifyEqual(length(get(gca, 'Children')), 10);
    
    tc.verifyEqual(length(findobj(gca, 'Type', 'line')), 5);
    tc.verifyEqual(length(findobj(gca, 'Type', 'text')), 5);
    labels = findobj(gca, 'Type', 'text');
    for i=1:5
        tc.verifyEqual(labels(i).Color, [0 1 0]);
        tc.verifyEqual(labels(i).FontWeight, 'normal');
    end
    
    % specify font color
    clf; plot_point(points, 'sequence', 'textcolor', [0.2 0.3 0.4]);
    tc.verifyEqual(length(get(gca, 'Children')), 10);
    
    tc.verifyEqual(length(findobj(gca, 'Type', 'line')), 5);
    tc.verifyEqual(length(findobj(gca, 'Type', 'text')), 5);
    labels = findobj(gca, 'Type', 'text');
    for i=1:5
        tc.verifyEqual(labels(i).Color, [0.2 0.3 0.4]);
        tc.verifyEqual(labels(i).FontWeight, 'normal');
    end
    
    % specify font weight
    clf; plot_point(points, 'sequence', 'bold');
    tc.verifyEqual(length(get(gca, 'Children')), 10);
    
    tc.verifyEqual(length(findobj(gca, 'Type', 'line')), 5);
    tc.verifyEqual(length(findobj(gca, 'Type', 'text')), 5);
    labels = findobj(gca, 'Type', 'text');
    for i=1:5
        tc.verifyEqual(labels(i).FontWeight, 'bold');
    end
    
    % specify font size, weight, color
    clf; plot_point(points, 'sequence', 'textsize', 23, 'bold', 'textcolor', [0.2 0.3 0.4]);
    tc.verifyEqual(length(get(gca, 'Children')), 10);
    
    tc.verifyEqual(length(findobj(gca, 'Type', 'line')), 5);
    tc.verifyEqual(length(findobj(gca, 'Type', 'text')), 5);
    labels = findobj(gca, 'Type', 'text');
    for i=1:5
        tc.verifyEqual(labels(i).FontSize, 23);
        tc.verifyEqual(labels(i).FontWeight, 'bold');
        tc.verifyEqual(labels(i).Color, [0.2 0.3 0.4]);
    end
end

function plotpoly_test(tc)

    opt.animate = false;
    % 2d and 3D

    opt.tag = [];
    opt.axis = [];
    
    P = [1 3 5; 1 5 1];
    
    clf; plot_poly(P);
    tc.verifyEqual(length(get(gca, 'Children')), 1);
    lines = findobj(gca, 'Type', 'line');
    tc.verifyEqual(length(lines), 1);
    tc.verifyEqual(length(lines(1).XData), 4);
    tc.verifyEqual(length(lines(1).YData), 4);
    tc.verifyEqual(lines(1).Marker, 'none');
    tc.verifyEqual(lines(1).LineStyle, '-');
        
    clf; plot_poly(P, 'g--');
    tc.verifyEqual(length(get(gca, 'Children')), 1);
    lines = findobj(gca, 'Type', 'line');
    tc.verifyEqual(length(lines), 1);
    tc.verifyEqual(length(lines(1).XData), 4);
    tc.verifyEqual(length(lines(1).YData), 4);
    tc.verifyEqual(lines(1).Marker, 'none');
    tc.verifyEqual(lines(1).LineStyle, '--');
    tc.verifyEqual(lines(1).Color, [0 1 0]);
    
    clf; plot_poly(P, 'edgecolor', 'r');
    tc.verifyEqual(length(get(gca, 'Children')), 1);
    lines = findobj(gca, 'Type', 'line');
    tc.verifyEqual(length(lines), 1);
    tc.verifyEqual(length(lines(1).XData), 4);
    tc.verifyEqual(length(lines(1).YData), 4);
    tc.verifyEqual(lines(1).Marker, 'none');
    tc.verifyEqual(lines(1).LineStyle, '-');
    %tc.verifyEqual(lines(1).Color, [1 0 0]);
    
    clf; plot_poly(P, 'edgecolor', [0.2 0.3 0.4]);
    tc.verifyEqual(length(get(gca, 'Children')), 1);
    lines = findobj(gca, 'Type', 'line');
    tc.verifyEqual(length(lines), 1);
    tc.verifyEqual(length(lines(1).XData), 4);
    tc.verifyEqual(length(lines(1).YData), 4);
    tc.verifyEqual(lines(1).Marker, 'none');
    tc.verifyEqual(lines(1).LineStyle, '-');
    %tc.verifyEqual(lines(1).Color, [0.2 0.3 0.4]);

    clf; plot_poly(P, 'tag', 'bob');
    lines = findobj(gca, 'Type', 'line');
    %tc.verifyEqual(lines(1).Tag, 'bob');  % no tag for line type
    
    %-------- patch mode
    
    clf; plot_poly(P, 'fillcolor', 'r');
    tc.verifyEqual(length(get(gca, 'Children')), 1);
    patch = findobj(gca, 'Type', 'patch');
    tc.verifyEqual(length(patch), 1);
    
    tc.verifyEqual(patch.FaceColor, [1 0 0]);
    tc.verifyEqual(patch.FaceAlpha, 1);
    tc.verifyEqual(patch.EdgeColor, [0 0 0]);
    tc.verifyEqual(patch.LineStyle, '-');
    tc.verifySize(patch.Faces, [1 3]);
    tc.verifySize(patch.Vertices, [3 2]);

    clf; plot_poly(P, 'fillcolor', 'g', 'alpha', 0.5);
    tc.verifyEqual(length(get(gca, 'Children')), 1);
    patch = findobj(gca, 'Type', 'patch');
    tc.verifyEqual(length(patch), 1);
    
    tc.verifyEqual(patch.FaceColor, [0 1 0]);
    tc.verifyEqual(patch.FaceAlpha, 0.5);
    tc.verifyEqual(patch.EdgeColor, [0 0 0]);
    tc.verifyEqual(patch.LineStyle, '-');
    tc.verifySize(patch.Faces, [1 3]);
    tc.verifySize(patch.Vertices, [3 2]);

    clf; plot_poly(P, 'fillcolor', 'b', 'tag', 'bob');
    patch = findobj(gca, 'Type', 'patch');
    %tc.verifyEqual(patch.Tag, 'bob');  % no tag 
end

function plot_sphere_test(tc)
end

function plot_box_test(tc)
end

function plot_arrow_test(tc)
end

function plot_homline_test(tc)
end

% plot_ellipse
function ellipse2d_test(tc)
    clf
    
    E = diag([9 4]);
    
    plot_ellipse(E)
    plot_ellipse(E, [4 0], 'r--');
    plot_ellipse(E, [0 4], 'edgecolor', 'g');
    plot_ellipse(E, [4 4], 'fillcolor', 'g');
    plot_ellipse(E, [0 8 0.5], 'edgecolor', 'r', 'fillcolor', 'c', 'alpha', 0.5, 'LineWidth', 3);
    plot_ellipse(E, [4 8], 'b--', 'LineWidth', 3);
    plot_ellipse(E, 'confidence', 0.5);
    
    h = plot_ellipse(E);
    plot_ellipse(E, [1 1], 'alter', h);
    tc.verifyError(@() plot_ellipse(E, 'alter', 7), 'SMTB:plot_ellipse:badarg');
end

function ellipse2d_animate_test(tc)
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
function ellipse3d_test(tc)
    clf
    
    E = diag([9 4 6]);
    
    plot_ellipse(E)

    
    clf
    plot_ellipse(E, 'edgecolor', 'g');

    
    clf
    plot_ellipse(E, 'fillcolor', 'g');
    
    clf
    plot_ellipse(E, 'fillcolor', 'g', 'shadow');

    
    clf
    plot_ellipse(E, 'fillcolor', 'g', 'edgecolor', 'r', 'LineWidth', 2);
    

    plot_ellipse(E, [0 8 1], 'edgecolor', 'r', 'fillcolor', 'c');
    plot_ellipse(E, [4 8 1], 'LineWidth', 3);

    axis equal
end

function ellipse3d_animate_test(tc)
    clf
    axis([-4 4 -4 4 -4 4]);
    
    E = diag([9 4 6]);
    
    h = plot_ellipse(E, 'g');
    for x = circle([0 0], 1)
        plot_ellipse(E, [x; 0], 'alter', h)
        pause(0.1)
    end
end
